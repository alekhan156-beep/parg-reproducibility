from pathlib import Path
import os, sys, subprocess, shutil, json, importlib.util
p=Path('/home/aliwhoami/parg'); run=Path((p/'runs/LATEST').read_text())
os.environ['PATH']='/home/aliwhoami/miniconda3/envs/parg/bin:'+os.environ['PATH']
smoke=run/'validation'; smoke.mkdir(exist_ok=True)
for name in ('data','models','results'): (smoke/name).mkdir(exist_ok=True)
shutil.copy2(p/'models/lgbm_final.pkl',smoke/'models/lgbm_final.pkl')
import pickle
from rdkit.Chem import Descriptors
with open(smoke/'models/lgbm_final.pkl','rb') as f: blob=pickle.load(f)
blob['descriptor_names']=[n for n,f in Descriptors._descList]
with open(smoke/'models/lgbm_final.pkl','wb') as f: pickle.dump(blob,f)
with open(p/'data/zinc_instock.smi') as f: lines=[next(f) for i in range(24)]
(smoke/'data/zinc_instock.smi').write_text(''.join(lines)+'invalid_smiles TEST_FAILURE\n')
r=subprocess.run([sys.executable,str(run/'run_parg_project.py'),'--stage','screen'],cwd=smoke,capture_output=True,text=True)
assert r.returncode==0,r.stderr
import pandas as pd,numpy as np
sc=pd.read_csv(smoke/'results/screen_all.csv'); errors=pd.read_csv(smoke/'results/screen_failures.csv')
assert len(sc)==24 and len(errors)==1
spec=importlib.util.spec_from_file_location('runner',run/'run_parg_project.py'); runner=importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
X=np.vstack([runner.feature(s) for s in sc.smiles]); expected=blob['model'].predict(X)
assert np.allclose(expected,sc.predicted_pchembl)
os.chdir(p)
ns={}; exec(compile((p/'scripts/40_allposes.py').read_text().split('\ntags =')[0],'rmsd','exec'),ns)
checks=[]
tags=(p/'results/xdock_targets.txt').read_text().split()
for tag,rec,src,ref,code in [('6OAK_M0M',r,p/f'results/redock_{r}.pdbqt',p/'structures/6OAK_ligand.pdb','M0M') for r in ('dry','wet')]+[(t,r,p/f'results/xdock/{t}_{r}.pdbqt',p/f'structures/xdock/{t}.pdb',t.split('_')[1]) for t in tags for r in ('dry','wet','B')]:
    dst=smoke/f'{tag}_{rec}.pdb'
    r=subprocess.run(['obabel','-ipdbqt',str(src),'-opdb','-O',str(dst),'-f','1','-l','1'],capture_output=True)
    assert r.returncode==0
    tmpl=ns['template'](code); q=ns['generic_query'](tmpl); refmol=ns['load_pdb'](str(ref)); pose=ns['load_pdb'](str(dst))
    assert pose is not None and refmol.GetNumAtoms()==pose.GetNumAtoms()==tmpl.GetNumAtoms(),tag
    mr=refmol.GetSubstructMatches(q,uniquify=False,maxMatches=301); mp=pose.GetSubstructMatches(q,uniquify=False,maxMatches=301)
    assert 0<len(mr)<=300 and 0<len(mp)<=300,(tag,rec,len(mr),len(mp))
    checks.append({'tag':tag,'rec':rec,'rmsd':ns['rmsd'](refmol,pose,q),'reference_matches':len(mr),'pose_matches':len(mp)})
result={'screen_scored':len(sc),'screen_failure_captured':len(errors),'screen_prediction_check':'passed','mapping_checks':checks}
(smoke/'checks.json').write_text(json.dumps(result,indent=2))
print('PASS: screening 24 valid + 1 invalid; mapping preflight',len(checks),'poses')
