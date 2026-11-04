from pathlib import Path
import shutil, subprocess, concurrent.futures, time, json, argparse, tempfile
parser=argparse.ArgumentParser(description='Run lessons and exercises in isolated folders. Requires R/JAGS and packages from 00_setup.R.')
parser.add_argument('--output', default=str(Path(tempfile.gettempdir())/'BSUecostats-validation'))
parser.add_argument('--workers',type=int,default=2)
parser.add_argument('--cases',nargs='*',help='Optional case names; omit to run everything.')
args=parser.parse_args()
repo=Path(__file__).resolve().parents[1]; root=Path(args.output).resolve();root.mkdir(parents=True,exist_ok=True)
runner=Path(__file__).with_name('run_lesson.R').resolve()
files={p.name.split('_')[0]:p.name for p in repo.glob('0*.R')}
cases=[]
def add(name,key,patches=()):cases.append((name,files[key],patches))
for key in sorted(files):add(key+'_default',key)
add('01_small','01',[('M <- 50','M <- 10'),('R <- 8','R <- 2')])
add('01_prior','01',[('M <- 50','M <- 10'),('R <- 8','R <- 2'),('print(rbind(prior.table(50, 8), prior.table(500, 80)), row.names = FALSE)','print(prior.table(10, 2), row.names = FALSE)')])
add('01_catch','01',[('C <- 600','C <- 1200')])
add('02_small','02',[('M <- 50','M <- 15')])
add('02_lowp','02',[('p <- 0.15','p <- 0.03')])
add('03_small','03',[('n.fish <- 400','n.fish <- 50')])
add('03_zero','03',[('b1 <- -0.8','b1 <- 0')])
add('04_outage','04',[('trap.out <- 41:45','trap.out <- 31:55')])
add('04_small','04',[('M <- rep(30,','M <- rep(10,')])
add('04_outage_long','04',[('trap.out <- 41:45','trap.out <- 31:55'),('ni <- 60000; nb <- 30000','ni <- 240000; nb <- 120000')])
add('04_small_long','04',[('M <- rep(30,','M <- rep(10,'),('ni <- 60000; nb <- 30000','ni <- 480000; nb <- 240000')])
add('05_facility','05',[('n.facility <- 1500','n.facility <- 5000')])
add('05_small','05',[('n.tagged <- 6000','n.tagged <- 600')])
# Only change the prior in the no-GOJ model.
s=(repo/files['05']).read_text(); at=s.index('# 8.'); tail=s[at:]; assert 'psi ~ dbeta(1, 1)' in tail
add('05_prior','05',[(tail,tail.replace('psi ~ dbeta(1, 1)','psi ~ dbeta(1, 9)'))])
add('06_more','06',[('log(8), log(150)','log(8), log(1500)')])
add('06_variable','06',[('sigma_psi <- 0.3','sigma_psi <- 0.8')])
add('06_free','06',[('restrict.spill <- TRUE','restrict.spill <- FALSE')])
add('06_positive_restricted','06',[('beta <- -0.8','beta <- 0.4')])
add('06_positive_free','06',[('beta <- -0.8','beta <- 0.4'),('restrict.spill <- TRUE','restrict.spill <- FALSE')])
add('06_repeats','06',[('run.repeats <- FALSE','run.repeats <- TRUE')])
add('06_variable_long','06',[('sigma_psi <- 0.3','sigma_psi <- 0.8')])
add('07_equal','07',[('a.true <- c(wild = 0.18','a.true <- c(wild = 0.09')])
add('07_zero','07',[('night_eff <- 3.0','night_eff <- 0')])
add('07_small','07',[('n.tagged <- c(wild = 80','n.tagged <- c(wild = 20')])
add('07_repeats','07',[('run.repeats <- FALSE','run.repeats <- TRUE')])
add('07_small_long','07',[('n.tagged <- c(wild = 80','n.tagged <- c(wild = 20')])
add('07a_equal','07a',[('result <- run_borrow_simulation()','wk <- borrow_inputs(); wk$aW <- wk$aH\n  result <- run_borrow_simulation(wk)')])
add('07a_small','07a',[('result <- run_borrow_simulation()','wk <- borrow_inputs(); wk$NW <- pmax(1, round(wk$NW / 10))\n  result <- run_borrow_simulation(wk)')])
add('07a_departure','07a',[('result <- run_borrow_simulation()','result <- run_borrow_simulation(wild.ge.multiplier = 0.8)')])
add('08_small','08',[('n.pit <- 600','n.pit <- 100')])
add('08_fewer','08',[('\nB <- 2000','\nB <- 200')])
add('09_small','09',[('n.tagged <- 3000','n.tagged <- 500')])
add('09_more','09',[('n.reps <- 50','n.reps <- 200')])
if args.workers < 1: parser.error('--workers must be positive')
if args.cases:
 unknown=set(args.cases)-{c[0] for c in cases}
 if unknown: parser.error('Unknown cases: '+', '.join(sorted(unknown)))
def run(case):
 name,filename,patches=case; d=root/name;d.mkdir(exist_ok=True)
 # A failed rerun must never leave an old result looking current.
 for artifact in ("result.rds", "direct_fits.rds"):
  (d/artifact).unlink(missing_ok=True)
 for p in repo.glob('*.R'):shutil.copy2(p,d/p.name)
 s=(d/filename).read_text()
 for old,new in patches:
  if old not in s:raise ValueError(name+': missing '+old[:100])
  s=s.replace(old,new,1)
 if name == '07_small_long':
  s=s.replace('n.iter = 6000, n.burnin = 2000', 'n.iter = 24000, n.burnin = 8000')
 if name == '06_variable_long':
  s=s.replace('n.iter = 10000, n.burnin = 5000', 'n.iter = 40000, n.burnin = 20000')
 (d/filename).write_text(s)
 t=time.time();print('START '+name,flush=True)
 with (d/'run.log').open('w') as f:r=subprocess.run(['Rscript',str(runner),str(d),filename],stdout=f,stderr=subprocess.STDOUT)
 print(f'DONE {name} exit={r.returncode} seconds={time.time()-t:.1f}',flush=True)
 return dict(case=name,exit=r.returncode,seconds=round(time.time()-t,1))
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
 futures={ex.submit(run,c):c[0] for c in cases if args.cases is None or c[0] in args.cases};results=[]
 for f in concurrent.futures.as_completed(futures):
  try:results.append(f.result())
  except Exception as e:
   print('SCHEDULER ERROR '+str(e),flush=True)
   results.append(dict(case=futures[f],exit=1,error=str(e)))
  (root/'manifest.json').write_text(json.dumps(results,indent=2))

raise SystemExit(int(any(r['exit'] != 0 for r in results)))
