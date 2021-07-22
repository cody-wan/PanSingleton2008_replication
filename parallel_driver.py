from datetime import datetime
from operator import mul
import pathlib
from scipy.optimize import minimize
import main
import sys 
import multiprocessing
import pickle
import time

def func(country_i, method='trust-constr'):
    XGUESS = [0, -0.75, 2.5, 1.5, 5, 0.025, 0.025]
    
    XUB = [7, 0.5, 5, 10, 10, 0.05, 0.05]
    XLB = [-7, -2, 0.0001, -7, 0.0001, 0.00001, 0.00001]

    ###############################################################
    ################# solve parameters via MLE ##################
    ###############################################################

    # method = 'L-BFGS-B' # stuck in local min, slow
    # method = 'TNC' # slow
    # method = 'SLSQP' # unstable
    # method = 'trust-constr' # pretty close, but with default xtol, still stuck in local min

    method_options = {'SLSQP':dict(maxiter=5000, disp=True), 
                      'trust-constr': dict(maxiter=5000, verbose=2, xtol=1e-8, gtol=1e-8)}

    now = datetime.now()

    # logs/method/country/time_at_logging/[console, csv, res object]
    log_file_path = f"logs/{method}/{main.country_names[country_i]}/{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}"
    # make folder
    pathlib.Path(f"{log_file_path}").mkdir(parents=True, exist_ok=True)
    log_csv_path = f"{log_file_path}/log.csv"

    # redirect console output to local file (for record-keeping)
    sys.stdout = open(f"{log_file_path}/console.txt", "w")
    print(f"writing to {log_file_path}")
    

    # document start time
    with open(log_csv_path, "a") as f:
        f.write(f"\tstart time: {now}\n")

    res = minimize(main.FUNC, args=(country_i, True, log_csv_path, False, False), x0=XGUESS, method=method, bounds=main.Bounds(lb=XLB, ub=XUB), options=method_options[method])

    with open(log_csv_path, "a") as f:
        f.write(f"\n\tend time: {datetime.now()}\n")

    sys.stdout.close()

    # save res object
    with open(f"{log_file_path}/res.pkl", "wb") as f:
        pickle.dump(res, f, pickle.HIGHEST_PROTOCOL)

# def func(i, m):
#     print(i)
#     time.sleep(5)

if __name__ == "__main__":

    # run concurrent func
    methods = ['trust-constr', 'SLSQP']
    # methods = ['trust-constr']
    country_nums = [i for i in range(15)]

    processes = [multiprocessing.Process(target=func, args=(country_i, method,)) for method in methods for country_i in country_nums]
    for p in processes:
        p.start()
    for p in processes:
        p.join()




