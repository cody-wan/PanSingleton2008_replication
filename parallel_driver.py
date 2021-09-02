from datetime import datetime
from operator import mul
import pathlib
from scipy.optimize import minimize
import main
import sys 
import numpy as np
import pandas as pd 
import multiprocessing
import pickle
import time
import resource


XGUESSES = {'Brazil': [3.52, -0.83, 0.75, -3.38, 0.69, 0.0059, 0.0035],
            'Bulgaria': [-1.15, 0.18, 1.58, -1.97, 0.34, 0.0017, 0.0013],
            'Chile': [-0.43, 0.06, 1.41, -2.82, 0.41, 0.0006, 0.0008],
            'China': [-0.50, 0.07, 1.12, -1.83, 0.29, 0.0004, 0.0005],
            'Colombia': [4.50, -1.16, 0.34, -5.54, 1.32, 0.0055, 0.0051],
            'Croatia': [-0.94, 0.14, 1.43, -3.37, 0.56, 0.0011, 0.0007],
            'Korea': [-0.50, 0.07, 1.09, -4.88, 0.81, 0.0007, 0.0006],
            'Malaysia': [-0.72, 0.13, 1.10, -1.07, 0.18, 0.0008, 0.0008],
            'Mexico': [-0.65, 0.11, 1.44, -3.86, 0.65, 0.0016, 0.0011],
            'Philippines': [0.17, -0.03, 1.05, -1.45, 0.31, 0.0054, 0.0037],
            'Poland': [-0.38, 0.05, 1.15, -4.47, 0.69, 0.0004, 0.0003],
            'Russia': [-1.78, 0.25, 2.20, -1.86, 0.15, 0.0024, 0.0023],
            'SAfrica': [-0.47, 0.06, 1.65, -2.82, 0.45, 0.0012, 0.0010],
            'Thailand': [-0.47, 0.08, 0.96, -4.39, 0.78, 0.0009, 0.0009],
            'Turkey': [-0.98, 0.18, 1.56, -1.50, 0.26, 0.0060, 0.0026]}



METHOD_OPTIONS = {'SLSQP':dict(maxiter=10000, disp=True, ftol=1e-12), 
                  'trust-constr': dict(maxiter=10000, verbose=2, xtol=1e-12, gtol=1e-12),
                  'L-BFGS-B': dict(maxfun=15000, maxiter=15000, iprint=101, ftol=1e-12, gtol=1e-6, eps=1e-10),
                  'TNC': dict(maxfun=5000)
                  }


def func(CDS1Y, CDS3Y, CDS5Y, PFILE, country_name, BATCH_NUM, XB_FLAG, LN_FLAG, method='trust-constr'):
    # if country in AEJ, use AEJ result as starting guess
    # if not, use average of AEJ results
    XGUESS = np.array(XGUESSES.get(country_name, np.array([-0.052, -0.043,  1.26, -3.014, 0.53,  0.0023,  0.0017])))
    
    XUB = [7, 0.5, 5, 10, 10, 0.05, 0.05]
    XLB = [-7, -2, 0.0001, -7, 0.0001, 0.00001, 0.00001]

    ###############################################################
    ################# solve parameters via MLE ##################
    ###############################################################

    # method = 'L-BFGS-B' # stuck in local min, slow
    # method = 'TNC' # slow
    # method = 'SLSQP' # unstable
    # method = 'trust-constr' # pretty close, but with default xtol, still stuck in local min

    now = datetime.now()

    log_file_path = f"logs/{DATA_V}/{BATCH_NUM}/{method}/{country_name}/{now.year}_{now.month}_{now.day}"
    # make folder
    pathlib.Path(log_file_path).mkdir(parents=True, exist_ok=True)
    log_csv_path = f"{log_file_path}/log.csv"

    # redirect console output to local file (for record-keeping)
    sys.stdout = open(f"{log_file_path}/console.txt", "w")
    print(f"writing to {log_file_path}")
    

    with open(log_csv_path, "a") as f:
        f.write(f"\tstart time: {now}\n")

    res = minimize(main.FUNC_NLLK, args=(CDS1Y, CDS3Y, CDS5Y, PFILE, country_name, XB_FLAG, LN_FLAG, True, log_csv_path, False, False), 
                              x0=XGUESS, method=method, bounds=main.Bounds(lb=XLB, ub=XUB), 
                              options=METHOD_OPTIONS[method])

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

    # read specifications:
    methods = ['trust-constr', 'SLSQP']
    # country_nums = [i for i in range(15)]

    batch_file = sys.argv[1]
    DATA_V = sys.argv[2]
    country_nums = [i for i in range(int(sys.argv[3]))] # 41 MKT; 15 AEJ

    BATCH_NUM, XB_FLAG, LN_FLAG = pd.read_csv(f"inputs/{DATA_V}/{batch_file}").iloc[0]

    # exit(0)
    processes = list()
    for country_i in country_nums:
        for method in methods:
            CDS1Y, CDS3Y, CDS5Y, PFILE, country_name = main.get_input(country_i, DATA_V=DATA_V)
            if len(PFILE) < 50:
                continue
            processes.append(multiprocessing.Process(target=func, args=(CDS1Y, CDS3Y, CDS5Y, PFILE, country_name, BATCH_NUM, XB_FLAG, LN_FLAG, method,)))

    for p in processes:
        p.start()
    for p in processes:
        p.join()




