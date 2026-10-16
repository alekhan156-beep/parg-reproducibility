"""Launch with the PARG conda Python inside WSL. Outputs are isolated per run."""
from pathlib import Path
import os, sys, json, shutil, subprocess, hashlib, datetime, time, traceback

SOURCE = Path('/home/aliwhoami/parg')
ENV = Path('/home/aliwhoami/miniconda3/envs/parg/bin')
os.environ['PATH'] = str(ENV) + ':' + os.environ.get('PATH', '')
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '4'
os.environ['MPLBACKEND'] = 'Agg'

def digest(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def save_json(p, obj):
    p = Path(p); tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, default=str)); tmp.replace(p)

def setup():
    run = SOURCE/'runs'/datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    run.mkdir(parents=True)
    for name in ('scripts', 'structures'):
        shutil.copytree(SOURCE/name, run/name)
    for name in ('data', 'models', 'results', 'figures', 'logs', 'configs'):
        (run/name).mkdir()
    for name in ('parg_raw.csv', 'fixed_split.json', 'chemprop_input.csv', 'zinc_instock.smi'):
        shutil.copy2(SOURCE/'data'/name, run/'data'/name)
    shutil.copytree(SOURCE/'data/ccd', run/'data/ccd')
    shutil.copy2(SOURCE/'results/xdock_targets.txt', run/'results/xdock_targets.txt')
    shutil.copy2(SOURCE/'vina_1.2.5_linux_x86_64', run/'vina_1.2.5_linux_x86_64')
    for seed in range(3):
        shutil.copy2(SOURCE/f'models/chemprop_seed{seed}/config.toml', run/f'configs/chemprop_seed{seed}.toml')
    shutil.copy2(__file__, run/'run_parg_project.py')
    save_json(run/'input_manifest.json', {str(p.relative_to(run)): {'sha256':digest(p),'bytes':p.stat().st_size}
              for parent in ('data','structures','scripts','configs') for p in (run/parent).rglob('*') if p.is_file()})
    save_json(run/'run_plan.json', {'source':str(SOURCE), 'inputs':'Frozen raw ChEMBL and ZINC snapshots; existing prepared receptors and ligand inputs copied and hashed.',
        'fresh_outputs':True, 'steps':['curation','splits','ML benchmarks and final model','Chemprop seeds 0,1,2 (50 epochs each)', '3.2M library screening','novelty','native redocking','dry/wet crossdocking','state B and box16 diagnostics','top100/random100 control','statistics and figures'],
        'limitations':['No new ChEMBL/ZINC download; reproduce saved dataset.', 'Prepared receptor chemistry is reused, not independently revalidated by this execution.', 'RMSD uses saved CCD connectivity query without pose alignment; bond orders and stereochemistry are not enforced.', 'No Boltz or molecular dynamics workflow exists in this project.']})
    with open(run/'launcher.log','w') as log:
        proc = subprocess.Popen([sys.executable,'-u',str(run/'run_parg_project.py'),'--worker'],cwd=run,
                                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    save_json(run/'launch.json',{'pid':proc.pid,'run':str(run),'started':datetime.datetime.now().isoformat()})
    (SOURCE/'runs/LATEST').write_text(str(run))
    print(json.dumps({'pid':proc.pid,'run':str(run)},indent=2))

def execute(name, cmd, timeout=86400):
    status['current_stage'] = name; status['updated'] = datetime.datetime.now().isoformat()
    save_json('status.json',status)
    print('START',name,flush=True); start=time.time()
    with open(Path('logs')/(name+'.log'),'w') as out:
        p=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,timeout=timeout)
    row={'stage':name,'returncode':p.returncode,'seconds':round(time.time()-start,2)}
    status['stages'].append(row); save_json('status.json',status)
    print('FINISH',row,flush=True)
    if p.returncode: raise RuntimeError(f'{name} failed; see logs/{name}.log')

def ml():
    import numpy as np, pandas as pd, pickle
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    from lightgbm import LGBMRegressor
    from sklearn.metrics import mean_squared_error,r2_score
    RDLogger.DisableLog('rdApp.*')
    df=pd.read_csv('data/parg_split.csv'); original=pd.read_csv(SOURCE/'data/parg_split.csv')
    assert df[['smiles','label','split_scaffold']].equals(original[['smiles','label','split_scaffold']]), 'Curated data/split differs from saved run'
    cp=pd.read_csv('data/chemprop_input.csv'); assert cp.smiles.tolist()==df.smiles.tolist()
    assert np.allclose(cp.label,df.label)
    splits=json.loads(Path('data/fixed_split.json').read_text())[0]
    sets={k:set(v) for k,v in splits.items()}
    assert sets['test']==set(np.flatnonzero(df.split_scaffold.eq('test')))
    assert not (sets['train'] & sets['val'] or sets['train'] & sets['test'] or sets['val'] & sets['test'])
    assert set.union(*sets.values())==set(range(len(df)))
    assert not (set(df.loc[df.split_scaffold.eq('train'),'scaffold']) & set(df.loc[df.split_scaffold.eq('test'),'scaffold']))
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    mols=[Chem.MolFromSmiles(s) for s in df.smiles]
    fps=[gen.GetFingerprint(m) for m in mols]
    desc=[n for n,f in Descriptors._descList]
    X=np.nan_to_num(np.array([np.concatenate([np.array(f),[fn(m) for n,fn in Descriptors._descList]]) for m,f in zip(mols,fps)]),nan=0,posinf=0,neginf=0)
    y=df.label.to_numpy(); metrics=[]; predictions=[]
    for split in ('split_random','split_scaffold'):
        tr=np.flatnonzero(df[split].eq('train')); te=np.flatnonzero(df[split].eq('test'))
        baseline=[]
        for i in te:
            sims=np.array(DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in tr]))
            top=np.argsort(sims)[-5:]; w=sims[top]
            baseline.append(np.average(y[tr][top],weights=w) if w.sum() else y[tr].mean())
        outputs=[('training_mean',-1,np.repeat(y[tr].mean(),len(te))),('tanimoto_5nn',-1,baseline)]
        for seed in range(5):
            model=LGBMRegressor(n_estimators=500,learning_rate=.05,num_leaves=31,random_state=seed,verbose=-1,n_jobs=4).fit(X[tr],y[tr])
            outputs.append(('lightgbm',seed,model.predict(X[te])))
        for name,seed,p in outputs:
            metrics.append(dict(split=split,model=name,seed=seed,n_train=len(tr),n_test=len(te),rmse=float(np.sqrt(mean_squared_error(y[te],p))),r2=r2_score(y[te],p)))
            predictions.extend(dict(split=split,model=name,seed=seed,index=int(i),smiles=df.smiles.iloc[i],truth=float(y[i]),prediction=float(v)) for i,v in zip(te,p))
    pd.DataFrame(metrics).to_csv('results/ml_metrics.csv',index=False)
    pd.DataFrame(predictions).to_csv('results/ml_predictions.csv',index=False)
    model=LGBMRegressor(n_estimators=500,learning_rate=.05,num_leaves=31,random_state=0,verbose=-1,n_jobs=4).fit(X,y)
    with open('models/lgbm_final.pkl','wb') as f: pickle.dump({'model':model,'n_desc':len(desc),'descriptor_names':desc},f)
    save_json('results/curation_counts.json',{'n':len(df),'active_ge6':int((y>=6).sum()),'inactive_le5':int((y<=5).sum()),'between5and6':int(((y>5)&(y<6)).sum()),'chemprop_split_sizes':{k:len(v) for k,v in sets.items()}})
    print(pd.DataFrame(metrics).to_string(index=False),flush=True)

def evaluate_chemprop():
    import pandas as pd, numpy as np
    from sklearn.metrics import mean_squared_error,r2_score
    truth=pd.read_csv('data/parg_split.csv'); truth=truth[truth.split_scaffold.eq('test')][['smiles','label']].rename(columns={'label':'truth'})
    metrics=pd.read_csv('results/ml_metrics.csv'); extra=[]
    for seed in range(3):
        files=list(Path(f'models/chemprop_seed{seed}').rglob('test_predictions.csv')); assert len(files)==1
        pred=pd.read_csv(files[0]).rename(columns={'label':'prediction'})
        m=truth.merge(pred,on='smiles',validate='one_to_one'); assert len(m)==len(truth)
        extra.append(dict(split='split_scaffold',model='chemprop',seed=seed,n_train=len(json.loads(Path('data/fixed_split.json').read_text())[0]['train']),n_test=len(m),rmse=float(np.sqrt(mean_squared_error(m.truth,m.prediction))),r2=r2_score(m.truth,m.prediction)))
    pd.concat([metrics,pd.DataFrame(extra)]).to_csv('results/ml_metrics.csv',index=False)
    print(extra)

def feature(smi):
    import numpy as np
    from rdkit import Chem,RDLogger
    from rdkit.Chem import Descriptors,rdFingerprintGenerator
    RDLogger.DisableLog('rdApp.*')
    try:
        m=Chem.MolFromSmiles(smi)
        if m is None: return None
        fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048).GetFingerprint(m)
        return np.nan_to_num(np.concatenate([np.array(fp),[f(m) for n,f in Descriptors._descList]]),nan=0,posinf=0,neginf=0).astype(np.float32)
    except Exception: return None

def screen():
    import pickle,csv,itertools,heapq,numpy as np
    from multiprocessing import Pool
    from rdkit.Chem import Descriptors
    with open('models/lgbm_final.pkl','rb') as f: blob=pickle.load(f)
    assert blob['descriptor_names']==[n for n,f in Descriptors._descList]
    model=blob['model']; model.set_params(n_jobs=4)
    heap=[]; seen=scored=0; start=time.time()
    with open('data/zinc_instock.smi') as src, open('results/screen_all.csv','w',newline='') as out, open('results/screen_failures.csv','w',newline='') as err, Pool(min(8,os.cpu_count() or 1)) as pool:
        w=csv.writer(out); w.writerow(['input_index','zinc_id','smiles','predicted_pchembl'])
        ew=csv.writer(err); ew.writerow(['input_index','input','reason'])
        while True:
            lines=list(itertools.islice(src,4096))
            if not lines: break
            batch=[]
            for line in lines:
                idx=seen; seen+=1; fields=line.split()
                if len(fields)<2: ew.writerow([idx,line.strip(),'malformed']); continue
                batch.append((idx,fields[0],fields[1]))
            feats=pool.map(feature,[x[1] for x in batch],chunksize=64)
            good=[j for j,f in enumerate(feats) if f is not None]
            if good:
                values=model.predict(np.vstack([feats[j] for j in good]))
                for j,value in zip(good,values):
                    idx,smi,zid=batch[j]; value=float(value); scored+=1
                    w.writerow([idx,zid,smi,value]); item=(value,-idx,zid,smi)
                    if len(heap)<5000: heapq.heappush(heap,item)
                    elif item>heap[0]: heapq.heapreplace(heap,item)
            for j,f in enumerate(feats):
                if f is None: ew.writerow([batch[j][0],batch[j][1],'parse_or_featurization_failure'])
            out.flush(); err.flush()
            save_json('screen_progress.json',{'processed':seen,'scored':scored,'seconds':time.time()-start,'molecules_per_second':seen/(time.time()-start)})
            print(f'processed={seen:,} scored={scored:,} speed={seen/(time.time()-start):.1f}/s',flush=True)
    with open('results/screen_ranked.csv','w',newline='') as out, open('results/screen_ranked.smi','w') as smiout:
        w=csv.writer(out); w.writerow(['rank','zinc_id','smiles','predicted_pchembl'])
        for rank,(value,idx,zid,smi) in enumerate(sorted(heap,reverse=True),1):
            w.writerow([rank,zid,smi,value]); smiout.write(f'{smi} {zid}\n')

def docking():
    import csv,random
    import numpy as np
    # Import only definitions from the archived all-pose implementation.
    ns={}; text=Path('scripts/40_allposes.py').read_text(); exec(compile(text.split('\ntags =')[0],'rmsd_definitions','exec'),ns)
    rows=[]; jobs=[]; failures=[]
    def call(cmd, logfile, timeout):
        with open(logfile,'w') as out:
            try: return subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,timeout=timeout).returncode
            except subprocess.TimeoutExpired: return 124
    def dock(label,rec,lig,dst,size=24,center=(27.20,2.84,.17)):
        dst=Path(dst); dst.parent.mkdir(parents=True,exist_ok=True)
        cmd=['./vina_1.2.5_linux_x86_64','--receptor',rec,'--ligand',lig,'--out',str(dst),'--exhaustiveness','16','--num_modes','9','--seed','42','--cpu','4']
        for ax,c in zip('xyz',center): cmd += ['--center_'+ax,str(c),'--size_'+ax,str(size)]
        rc=call(cmd,Path('logs')/(label+'.log'),3600)
        aff=ns['affinities'](dst) if rc==0 and dst.exists() else []
        jobs.append({'job':label,'command':cmd,'returncode':rc,'n_modes':len(aff)})
        save_json('results/docking_jobs.json',jobs)
        if not aff:
            failures.append(label); print('FAILED',label,rc,flush=True)
        else: print('DOCKED',label,aff[0],flush=True)
        return aff
    def measure(tag,rec,dst,reference,code,aff):
        tmpl=ns['template'](code); q=ns['generic_query'](tmpl); ref=ns['load_pdb'](reference)
        base=Path('results/tmp_poses')/(tag+'_'+rec)
        rc=call(['obabel','-ipdbqt',str(dst),'-opdb','-O',str(base)+'.pdb','-m'],Path('logs')/(tag+'_'+rec+'_convert.log'),120)
        if rc: raise RuntimeError('Pose conversion failed '+str(dst))
        for rank,a in enumerate(aff,1):
            pose=ns['load_pdb'](str(base)+str(rank)+'.pdb')
            if pose is None or ref.GetNumAtoms()!=tmpl.GetNumAtoms() or pose.GetNumAtoms()!=tmpl.GetNumAtoms(): raise RuntimeError('Incomplete RMSD atom coverage '+tag)
            mr=ref.GetSubstructMatches(q,uniquify=False,maxMatches=301); mp=pose.GetSubstructMatches(q,uniquify=False,maxMatches=301)
            if not mr or not mp or len(mr)>300 or len(mp)>300: raise RuntimeError('Missing or capped atom mapping '+tag)
            val=ns['rmsd'](ref,pose,q)
            if val is None: raise RuntimeError('Unmapped RMSD '+tag)
            rows.append(dict(ligand=tag,receptor=rec,rank=rank,affinity=a,true_rmsd=val))
        with open('results/all_docking_poses.csv','w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['ligand','receptor','rank','affinity','true_rmsd']); w.writeheader(); w.writerows(rows)
    for rec in ('dry','wet'):
        dst=f'results/redock_{rec}.pdbqt'; aff=dock('redock_'+rec,f'structures/rec_{rec}.pdbqt','structures/lig.pdbqt',dst)
        if aff: measure('6OAK_M0M','redock_'+rec,dst,'structures/6OAK_ligand.pdb','M0M',aff)
    tags=Path('results/xdock_targets.txt').read_text().split(); assert len(tags)==17
    for tag in tags:
        for rec in ('dry','wet','B'):
            dst=f'results/xdock/{tag}_{rec}.pdbqt'; aff=dock(tag+'_'+rec,f'structures/rec_{rec}.pdbqt',f'structures/xdock/{tag}.pdbqt',dst)
            if aff: measure(tag,rec,dst,f'structures/xdock/{tag}.pdb',tag.split('_')[1],aff)
        ref=ns['load_pdb'](f'structures/xdock/{tag}.pdb'); center=np.round(ref.GetConformer().GetPositions().mean(0),2)
        for rec in ('wet','B'):
            dst=f'results/box16/{tag}_{rec}.pdbqt'; aff=dock(tag+'_box16_'+rec,f'structures/rec_{rec}.pdbqt',f'structures/xdock/{tag}.pdbqt',dst,16,center)
            if aff: measure(tag,'box16_'+rec,dst,f'structures/xdock/{tag}.pdb',tag.split('_')[1],aff)
    if failures: raise RuntimeError('Docking benchmark failures: '+str(failures))
    # Fresh ligand preparation; preserve every failure in the denominator.
    top=[tuple(l.split()[:2]) for l in Path('results/screen_ranked.smi').read_text().splitlines()][:100]
    ids={z for s,z in top}
    library=[tuple(l.split()[:2]) for l in open('data/zinc_instock.smi') if len(l.split())>=2 and l.split()[1] not in ids]
    rnd=random.Random(42).sample(library,100); del library
    control=[]; Path('ligands').mkdir(exist_ok=True)
    for i in range(100):
        for group,(smi,zid) in [('top',top[i]),('random',rnd[i])]:
            lig=f'ligands/{zid}.pdbqt'; dst=f'results/control/{zid}.pdbqt'
            rc=call(['obabel',f'-:{smi}','-opdbqt','-O',lig,'--gen3d','-h'],Path('logs')/(zid+'_prep.log'),600)
            a=[]; reason='preparation_failed'
            if rc==0 and Path(lig).exists() and Path(lig).stat().st_size:
                a=dock('control_'+zid,'structures/rec_wet.pdbqt',lig,dst,16); reason='ok' if a else 'docking_failed'
            control.append(dict(group=group,zinc_id=zid,smiles=smi,status=reason,prep_returncode=rc,vina_kcal=a[0] if a else ''))
            with open('results/control_all.csv','w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(control[0])); w.writeheader(); w.writerows(control)
    print('Control status', {s:sum(r['status']==s for r in control) for s in {r['status'] for r in control}},flush=True)

def report():
    import pandas as pd, numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import mannwhitneyu,binomtest
    metrics=pd.read_csv('results/ml_metrics.csv'); poses=pd.read_csv('results/all_docking_poses.csv'); control=pd.read_csv('results/control_all.csv')
    summary={}
    for rec,g in poses.groupby('receptor'):
        first=g[g['rank'].eq(1)]; anypose=g.groupby('ligand').true_rmsd.min()
        summary[rec]={'n':len(first),'top1_under2':int(first.true_rmsd.lt(2).sum()),'sampled_under2':int(anypose.lt(2).sum())}
    cross=poses[poses.receptor.isin(['dry','wet'])]
    cross.to_csv('results/crossdock_allposes.csv',index=False)
    paired=cross[cross['rank'].eq(1)].pivot(index='ligand',columns='receptor',values='true_rmsd')
    a=int(((paired.dry<2)&(paired.wet>=2)).sum()); b=int(((paired.wet<2)&(paired.dry>=2)).sum())
    summary['paired_top1_exact_p']=binomtest(a,a+b,.5).pvalue if a+b else 1.
    ok=control[control.status.eq('ok')]; t=ok.loc[ok.group.eq('top'),'vina_kcal']; r=ok.loc[ok.group.eq('random'),'vina_kcal']
    summary['control']={'attempts':len(control),'top_success':len(t),'random_success':len(r),'top_mean':t.mean(),'random_mean':r.mean(),'mannwhitney_two_sided':mannwhitneyu(t,r,alternative='two-sided').pvalue if len(t) and len(r) else None}
    save_json('results/summary.json',summary)
    fig,ax=plt.subplots(figsize=(8,4)); g=metrics[metrics.split.eq('split_scaffold')].groupby('model').rmse.agg(['mean','std']).fillna(0)
    ax.bar(g.index,g['mean'],yerr=g['std']); ax.set_ylabel('Scaffold test RMSE (pChEMBL)'); fig.tight_layout(); fig.savefig('figures/ml_benchmark.png',dpi=200); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5)); x=np.arange(len(paired)); ax.scatter(x,paired.dry,label='Dry'); ax.scatter(x,paired.wet,label='Wet'); ax.axhline(2,color='grey',ls='--'); ax.set_xticks(x,paired.index,rotation=90); ax.set_ylabel('Top-ranked pose RMSD (Å)'); ax.legend(); fig.tight_layout(); fig.savefig('figures/crossdock.png',dpi=200); plt.close(fig)
    Path('RUN_REPORT.md').write_text('# PARG fresh execution\n\nCompleted from frozen input snapshots. See input_manifest.json and run_plan.json for provenance and scope.\n\n'+metrics.to_csv(index=False)+'\n\nDocking summary:\n\n```json\n'+json.dumps(summary,indent=2)+'\n```\n\nControl failures remain in control_all.csv; complete-case statistics do not resolve missing-outcome bias. RMSD preserves the archived connectivity-only mapping method and does not enforce stereochemistry.\n')

def worker():
    global status
    status={'state':'running','pid':os.getpid(),'stages':[],'started':datetime.datetime.now().isoformat()}; save_json('status.json',status)
    try:
        execute('01_curate',[sys.executable,'-u','scripts/02_curate.py'])
        execute('02_split',[sys.executable,'-u','scripts/03_split.py'])
        execute('03_ml',[sys.executable,'-u',__file__,'--stage','ml'])
        for seed in range(3):
            execute(f'04_chemprop_seed{seed}',[str(ENV/'chemprop'),'train','--config-path',f'configs/chemprop_seed{seed}.toml','--output-dir',f'models/chemprop_seed{seed}'])
        execute('05_chemprop_evaluation',[sys.executable,'-u',__file__,'--stage','evaluate_chemprop'])
        execute('06_screen',[sys.executable,'-u',__file__,'--stage','screen'])
        execute('07_novelty',[sys.executable,'-u','scripts/33_novelty.py'])
        execute('08_docking',[sys.executable,'-u',__file__,'--stage','docking'],timeout=7*86400)
        execute('09_report',[sys.executable,'-u',__file__,'--stage','report'])
        status.update(state='complete',current_stage=None)
    except Exception:
        status.update(state='failed',error=traceback.format_exc()); print(status['error'],flush=True)
    status['updated']=datetime.datetime.now().isoformat(); save_json('status.json',status)

if __name__=='__main__':
    if '--worker' in sys.argv: worker()
    elif '--stage' in sys.argv: globals()[sys.argv[sys.argv.index('--stage')+1]]()
    else: setup()
