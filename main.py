from numba.np.ufunc import parallel
import numpy as np 
import pandas as pd 
import pathlib 
import numba
from scipy.linalg import solve_banded
from scipy.optimize import minimize, Bounds
from datetime import datetime
from scipy import sparse

# countries to be included in estimation:
TOTAL_COUNTRY_NAMES = ['Brazil', 'Bulgaria', 'Chile', 'China', 'Colombia', 
                       'Croatia', 'Korea', 'Malaysia', 'Mexico', 'Philippines', 
                       'Poland', 'Russia', 'SAfrica', 'Thailand', 'Turkey']


# helper functions
@numba.jit(nopython=True)
def diagonal_form(a, upper=1, lower=1):
    """
        a is a numpy square matrix
        this function converts a square matrix to diagonal ordered form
        returned matrix in ab shape which can be used directly for scipy.linalg.solve_banded
    """
    n = np.shape(a)[1]
    # assert(np.all(a.shape ==(n,n)))
    
    ab = np.zeros((2*n-1, n))
    
    for i in np.arange(n):
        ab[i,(n-1)-i:] = np.diag(a,(n-1)-i)
        
    for i in np.arange(n-1): 
        ab[(2*n-2)-i,:i+1] = np.diag(a,i-(n-1))

    mid_row_inx = np.int(ab.shape[0]/2)

    upper_rows = np.array([mid_row_inx - i for i in np.arange(start=1, stop=upper+1)])
    upper_rows = np.flip(upper_rows)
    upper_rows = np.append(upper_rows, mid_row_inx)

    lower_rows = np.array([mid_row_inx + i for i in np.arange(start=1, stop=lower+1)])

    keep_rows = np.append(upper_rows, lower_rows) 

    return ab[keep_rows,:]


@numba.jit(nopython=True)
def construct_B(ALP, BET, XSTEP, TSTEP, SIG, X, XB_FLAG):
    """
        XB_FLAG: T/F; T: use X[99] same as fortran code, F: use X[199], mathematically correct
    """

    B = np.zeros((200,200))
    # edge values
    B[0, 0] = -(ALP - BET * X[0]) / XSTEP - np.exp(X[0]) - 1/TSTEP
    B[0, 1] = (ALP - BET * X[0]) / XSTEP

    if XB_FLAG:
        B[199, 198] = -(ALP - BET * X[99]) / XSTEP
        B[199, 199] = (ALP - BET * X[99]) / XSTEP - np.exp(X[99]) - 1/TSTEP
    else:
        B[199, 198] = -(ALP - BET * X[199]) / XSTEP
        B[199, 199] = (ALP - BET * X[199]) / XSTEP - np.exp(X[199]) - 1/TSTEP


    for I in range(1, 199):
        B[I, I-1] = SIG**2 / (2 * XSTEP**2) - (ALP - BET*X[I]) / (2 * XSTEP)
        B[I, I] = -SIG**2 / XSTEP**2 - np.exp(X[I]) - 1/TSTEP
        B[I, I+1] = SIG**2 / (2 * XSTEP**2) + (ALP - BET*X[I]) / (2 * XSTEP)
    
    return B


def PDE(PARM, init_cond, measure, XB_FLAG):
    """
        PDE: solve PDE via finite difference
        args:
            measure: 'Q' or 'P'
    """
    XLOWER = 0.0001
    XUPPER = 0.25

    TSTEP = 1./12.
    XSTEP = (np.log(XUPPER) - np.log(XLOWER))/200

    if measure == 'Q':
        ALP = PARM[0]
        BET = PARM[1]
    elif measure == 'P':
        ALP = PARM[3] # ALPH
        BET = PARM[4] # BETH

    SIG = PARM[2]

    X = np.array([np.log(XLOWER) + XSTEP * float(i) for i in np.arange(200)])
    XX = np.exp(X)

    if init_cond == 0:
        # initial condition, F[i, 1] = 1 for i in 0:199. F[:,1+] are populated iteratively
        mat = np.ones((200,61)) 
    elif init_cond == 1:
        # initial condition, G[i, 1] = exp(A[i]) for i in 0:199. G[:,1+] are populated iteratively
        mat = np.ones((200,61)) 
        mat[:,0] = XX


    for IT in np.arange(start=1, stop=61):

        A = np.array([ -1 * mat[i, IT-1] / TSTEP for i in np.arange(200)])
        B = construct_B(ALP, BET, XSTEP, TSTEP, SIG, X, XB_FLAG)

        # solve for x in Bx=A
        # 1. solve via solve_banded()
        ab = diagonal_form(B)
        res = solve_banded((1,1), ab, A)

        # 2. solve via least square; takes twice as much time as #1
        # res = np.linalg.lstsq(B, A)[0]

        for i in np.arange(200):
            mat[i, IT] = res[i]

    return mat, XX

@numba.jit(nopython=True)
def construct_CDS(PFILE, WR, ICAR_2, F, G):

    CDS1T = np.zeros(200)
    CDS3T = np.zeros(200)
    CDS5T = np.zeros(200)

    for i in np.arange(200):

        CDS1_sum1, CDS1_sum2 = 0, 0
        CDS3_sum1, CDS3_sum2 = 0, 0
        CDS5_sum1, CDS5_sum2 = 0, 0

        for j in np.arange(1, 20+1):
            # 1 YR CDS
            if j <= 4:
                CDS1_sum1 += PFILE[ICAR_2, 13*j - 1] * F[i, 1+j*3 - 1]
                CDS1_sum2 += PFILE[ICAR_2, 13*j - 1] * G[i, 1+j*3 - 1]
            # 3 YR CDS
            if j <= 12:
                CDS3_sum1 += PFILE[ICAR_2, 13*j - 1] * F[i, 1+j*3 - 1]
                CDS3_sum2 += PFILE[ICAR_2, 13*j - 1] * G[i, 1+j*3 - 1]
            # 5 YR CDS
            if j <= 20:
                CDS5_sum1 += PFILE[ICAR_2, 13*j - 1] * F[i, 1+j*3 - 1]
                CDS5_sum2 += PFILE[ICAR_2, 13*j - 1] * G[i, 1+j*3 - 1]

        CDS1T[i] = WR * CDS1_sum2 / CDS1_sum1
        CDS3T[i] = WR * CDS3_sum2 / CDS3_sum1
        CDS5T[i] = WR * CDS5_sum2 / CDS5_sum1
    
    return CDS1T, CDS3T, CDS5T


def EVAL(PFILE, PARM, ICAR_2, measure, XB_FLAG, WR=0.75):
    """
        EVAL: evaluates model value of CDS1, CDS3 and CDS5
        args: PARM: parameters
              ICAR_2: time indicator
    """    

    F, X = PDE(PARM, init_cond=0, measure=measure, XB_FLAG=XB_FLAG)
    G, _ = PDE(PARM, init_cond=1, measure=measure, XB_FLAG=XB_FLAG)

    CDS1T, CDS3T, CDS5T = construct_CDS(PFILE, WR, ICAR_2, F, G)

    return X, CDS1T, CDS3T, CDS5T

@numba.jit(nopython=True)
def XINTP(CDS3, X, CDS1T, CDS3T, CDS5T):
    """
        interpolates CDS prices to match market data under Q measure
    """
    # default values
    xlam=X[199]
    cds1x=CDS1T[199]
    cds3x=CDS3T[199]
    cds5x=CDS5T[199]
    # not exactly sure what the default value is in their fortran 77 code
    deriv=(CDS3T[199] - CDS3T[198]) / (np.log(X[199]) - np.log(X[198]))

    for i in np.arange(199):
        if ((CDS3 >= CDS3T[i]) and (CDS3 <= CDS3T[i+1])):
            weight = (CDS3 - CDS3T[i]) / (CDS3T[i+1] - CDS3T[i])
            xlam = (1 - weight) * X[i] + weight * X[i+1]
            cds1x = (1 - weight) * CDS1T[i] + weight * CDS1T[i+1]
            cds3x = (1 - weight) * CDS3T[i] + weight * CDS3T[i+1]
            cds5x = (1 - weight) * CDS5T[i] + weight * CDS5T[i+1]

            deriv = (CDS3T[i+1] - CDS3T[i]) / (np.log(X[i+1]) - np.log(X[i])) # NOTE X is evenly spaced out

            break
    
    return xlam, cds1x, cds3x, cds5x, deriv

@numba.jit(nopython=True)
def XINTPP(X, XLAM, CDS1T, CDS3T, CDS5T):
    """
        interpolates CDS prices to match market data under P measure
    """

    # default values
    xlam=X[199]
    cds1x=CDS1T[199]
    cds3x=CDS3T[199]
    cds5x=CDS5T[199]

    for i in np.arange(199):
        if ((XLAM >= X[i]) and (XLAM <= X[i+1])):
            weight = (XLAM - X[i]) / (X[i+1] - X[i])
            cds1x = (1 - weight) * CDS1T[i] + weight * CDS1T[i+1]
            cds3x = (1 - weight) * CDS3T[i] + weight * CDS3T[i+1]
            cds5x = (1 - weight) * CDS5T[i] + weight * CDS5T[i+1]

            break
    
    return cds1x, cds3x, cds5x


def FUNC_NLLK(PARM, CDS1Y, CDS3Y, CDS5Y, PFILE, country_name, 
                        XB_FLAG, LN_FLAG, 
                        logging=True, log_file_path='logs/log.csv', print_flag=False, CDSLog_flag=False, CDSLog_path="outputs"):
    """
        FUNC: solves CDS spreads PDE, interpolates model spreads, returns negative log likelihood
        args:
            PARM: theta; 
            ICAR_1: country indicator
            LN_FLAG: boolean; T: include `ln(derive)` in log-likelihood, F: no
    """

    # prepare data for FUNC_NLLK
    # 1. date-aligned CDS 1Y, 3Y, 5Y data

    resQ, resP = [], []

    ALP = PARM[0]
    BET = PARM[1]
    SIG = PARM[2]
    ALPH = PARM[3]
    BETH = PARM[4]
    SIG1 = PARM[5]
    SIG2 = PARM[6]

    if(SIG <= 0):
        SIG=0.0001
    if(SIG1 <= 0):
        SIG=0.00001
    if(SIG2 <= 0):
        SIG=0.00001
    if(BETH <= 0):
        SIG=0.001

    PARM = np.array([ALP, BET, SIG, ALPH, BETH, SIG1, SIG2])

    TOTAL, TOTALL = 0, 0
 
    for ICAR_2 in range(len(CDS1Y)):
        # read market CDS data
        cds1 = CDS1Y[ICAR_2]
        cds3 = CDS3Y[ICAR_2]
        cds5 = CDS5Y[ICAR_2]

        # get model CDS values under Q measure
        X, CDS1T, CDS3T, CDS5T = EVAL(PFILE, PARM, ICAR_2, measure='Q', XB_FLAG=XB_FLAG)
        # interpolate to match market data under Q measure
        xlam, cds1x, cds3x, cds5x, deriv = XINTP(cds3, X, CDS1T, CDS3T, CDS5T)
        resQ.append([xlam, cds1x, cds3x, cds5x, deriv])

        TOTALL += (cds1 - cds1x)**2 + (cds3 - cds3x)**2 + (cds5 - cds5x)**2

        # get model CDS values under P measure
        X, CDS1TP, CDS3TP, CDS5TP = EVAL(PFILE, PARM, ICAR_2, measure='P', XB_FLAG=XB_FLAG)
        # interpolate to match market data under P measure
        cds1xp, cds3xp, cds5xp = XINTPP(X, xlam, CDS1TP, CDS3TP, CDS5TP)
        resP.append([cds1xp, cds3xp, cds5xp])

        if(ICAR_2==0):
            XLAMM = xlam
        
        # compute log-likelihood
        ER1 = cds1 - cds1x
        ER2 = cds5 - cds5x

        V = SIG * SIG / (2 * BETH) * (1 - np.exp(-2 * BETH / 12))
        XM = np.log(XLAMM) * np.exp(-1 * BETH /12.) + ALPH/BETH * (1 - np.exp(-1 * BETH / 12.))

        # log likelihood
        if LN_FLAG:
            TOTAL = TOTAL -0.5 * np.log(2 * np.pi * SIG1**2) \
                        - ER1 * ER1 / (2 * SIG1**2) \
                        - 0.5 * np.log(2 * np.pi * SIG2**2) \
                        - ER2 * ER2 / (2 * SIG2**2) \
                        -0.5 * np.log(2 * np.pi * V) \
                        -1 * (np.log(xlam) - XM)**2 / (2 * V) -1 * np.log(deriv)
        else:
            TOTAL = TOTAL -0.5 * np.log(2 * np.pi * SIG1**2) \
                        - ER1 * ER1 / (2 * SIG1**2) \
                        - 0.5 * np.log(2 * np.pi * SIG2**2) \
                        - ER2 * ER2 / (2 * SIG2**2) \
                        -0.5 * np.log(2 * np.pi * V) \
                        -1 * (np.log(xlam) - XM)**2 / (2 * V) # -1 * np.log(deriv)
        
        XLAMM = xlam

    FVALL = TOTALL / (3 * len(CDS1Y))
    FVALL = FVALL**0.5
    FVAL = -TOTAL

    # write model CDS spreads
    if CDSLog_flag:
        XLAM, CDS1X, CDS3X, CDS5X, DERIV = np.array(resQ).T
        CDS1XP, CDS3XP, CDS5XP = np.array(resP).T
        CDS_df = pd.DataFrame.from_dict(dict(XLAM=XLAM, CDS1Y=CDS1Y, CDS1X=CDS1X, CDS1XP=CDS1XP, 
                                    CDS3=CDS3Y, CDS3X=CDS3X, CDS3XP=CDS3XP, 
                                    CDS5=CDS5Y, CDS5X=CDS5X, CDS5XP=CDS5XP))
        CDS_df.to_csv(f"{CDSLog_path}/{country_name}_CDS.csv", index=False)

    # write parameter value (at each function call/iteration)
    if logging:
        with open(log_file_path, "a") as f:
            now = datetime.now()
            f.write(f"{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute},")
            f.write(",".join(map(str,PARM)))
            f.write(f", {FVAL}, {FVALL}\n")
    # print out parameter value in console
    if print_flag:
        print(f'{PARM},{FVAL}, {FVALL}')

    return FVAL # negative log likelihood (to be minimized)



def find_overlap(ICAR_1, DATA_V):
    """
        return largest subset of data where all CDS1Y, CDS3Y, CDS5Y and PFILE data are available
    """

    CDS1Y = pd.read_csv(f'inputs/{DATA_V}/CDS1Y.csv', index_col=0, parse_dates=[0]).dropna(axis=1, how='all')/10000
    CDS3Y = pd.read_csv(f'inputs/{DATA_V}/CDS3Y.csv', index_col=0, parse_dates=[0]).dropna(axis=1, how='all')/10000
    CDS5Y = pd.read_csv(f'inputs/{DATA_V}/CDS5Y.csv', index_col=0, parse_dates=[0]).dropna(axis=1, how='all')/10000

    # NOTE: make sure first len(TOTAL_COUNTRY_NAMES) countries in CDS files are the same as TOTAL_COUNTRY_NAMES
    CDS1Y.columns = TOTAL_COUNTRY_NAMES + list(CDS1Y.columns[len(TOTAL_COUNTRY_NAMES):])
    CDS3Y.columns = TOTAL_COUNTRY_NAMES + list(CDS3Y.columns[len(TOTAL_COUNTRY_NAMES):])
    CDS5Y.columns = TOTAL_COUNTRY_NAMES + list(CDS5Y.columns[len(TOTAL_COUNTRY_NAMES):])

    country_names = [val for val in CDS1Y.columns.tolist() if val in CDS3Y.columns.tolist() and val in CDS5Y.columns.tolist()]

    if ICAR_1 > len(country_names):
        return 0, 0, 0, 0, "NA"

    # read 1Y, 3Y, 5Y market CDS data
    # take subset where data of all tenors are available
    country = country_names[ICAR_1]

    CDS1Y_t0 = CDS1Y[country].first_valid_index()
    CDS3Y_t0 = CDS3Y[country].first_valid_index()
    CDS5Y_t0 = CDS5Y[country].first_valid_index()

    CDS1Y_t1 = CDS1Y[country].last_valid_index()
    CDS3Y_t1 = CDS3Y[country].last_valid_index()
    CDS5Y_t1 = CDS5Y[country].last_valid_index()

    # PFILE
    PFILE = pd.read_csv(f"inputs/{DATA_V}/PFILE.csv", index_col=0, parse_dates=[0])
    P_t0, P_t1 = PFILE.index[0], PFILE.index[len(PFILE)-1]

    t0, t1 = max(max(max(CDS1Y_t0, CDS3Y_t0), CDS5Y_t0), P_t0), min(min(min(CDS1Y_t1, CDS3Y_t1), CDS5Y_t1), P_t1)

    return CDS1Y, CDS3Y, CDS5Y, PFILE, country, t0, t1


def get_input(ICAR_1, DATA_V):
    """
        prepare input for FUNC_NLLK

        return:
            1Y, 3Y, 5Y CDS list
            PFILE: corresponding discounting values
            country_name
    """

    if DATA_V == "AEJ_data":
        # read data and add country index header

        CDS1Y = pd.read_csv(f'inputs/{DATA_V}/CDS1Y.TXT', delimiter = "\t", index_col=False, names=[i for i in range(len(TOTAL_COUNTRY_NAMES))], dtype=float)/10000
        CDS3Y = pd.read_csv(f'inputs/{DATA_V}/CDS3Y.TXT', delimiter = "\t", index_col=False, names=[i for i in range(len(TOTAL_COUNTRY_NAMES))], dtype=float)/10000
        CDS5Y = pd.read_csv(f'inputs/{DATA_V}/CDS5Y.TXT', delimiter = "\t", index_col=False, names=[i for i in range(len(TOTAL_COUNTRY_NAMES))], dtype=float)/10000
        # read TPFILE
        PFILE = pd.read_csv(f'inputs/{DATA_V}/TPFILE', delimiter = "     ", header=None, dtype=float, engine='python').values


        return CDS1Y.values[:,ICAR_1], CDS3Y.values[:,ICAR_1], CDS5Y.values[:,ICAR_1], PFILE, TOTAL_COUNTRY_NAMES[ICAR_1]
    
    else: # "AEJ_data", "Test_data", etc
        CDS1Y, CDS3Y, CDS5Y, PFILE, country, t0, t1 = find_overlap(ICAR_1=ICAR_1, DATA_V=DATA_V)

        print(f"{country} {t0.year}-{t0.month}-{t0.day} - {t1.year}-{t1.month}-{t1.day} # obs: {len(PFILE.loc[t0:t1])}")

        return CDS1Y[country].loc[t0: t1].values, \
               CDS3Y[country].loc[t0: t1].values, \
               CDS5Y[country].loc[t0: t1].values, \
               PFILE.loc[t0: t1].values, \
               country



if __name__ == '__main__':


    # variables to be estimated, guesses, lower bound, upper bound
    # 0. alpha 1. beta 2. sigma 3. alpha_hat 4. beta_hat 5. sig_err1 6. eig_err2

    XGUESSES = [[3.52, -0.83, 0.75, -3.38, 0.69, 0.0059, 0.0035],
                [-1.15, 0.18, 1.58, -1.97, 0.34, 0.0017, 0.0013],
                [-0.43, 0.06, 1.41, -2.82, 0.41, 0.0006, 0.0008],
                [-0.50, 0.07, 1.12, -1.83, 0.29, 0.0004, 0.0005],
                [4.50, -1.16, 0.34, -5.54, 1.32, 0.0055, 0.0051],
                [-0.94, 0.14, 1.43, -3.37, 0.56, 0.0011, 0.0007],
                [-0.50, 0.07, 1.09, -4.88, 0.81, 0.0007, 0.0006],
                [-0.72, 0.13, 1.10, -1.07, 0.18, 0.0008, 0.0008],
                [-0.65, 0.11, 1.44, -3.86, 0.65, 0.0016, 0.0011],
                [0.17, -0.03, 1.05, -1.45, 0.31, 0.0054, 0.0037],
                [-0.38, 0.05, 1.15, -4.47, 0.69, 0.0004, 0.0003],
                [-1.78, 0.25, 2.20, -1.86, 0.15, 0.0024, 0.0023],
                [-0.47, 0.06, 1.65, -2.82, 0.45, 0.0012, 0.0010],
                [-0.47, 0.08, 0.96, -4.39, 0.78, 0.0009, 0.0009],
                [-0.98, 0.18, 1.56, -1.50, 0.26, 0.0060, 0.0026]]

    XGUESS = [0, -0.75, 2.5, 1.5, 5, 0.025, 0.025]
    
    XUB = [7, 0.5, 5, 10, 10, 0.05, 0.05]
    XLB = [-7, -2, 0.0001, -7, 0.0001, 0.00001, 0.00001]

    for country_i in range(1):
        CDS1Y, CDS3Y, CDS5Y, PFILE, country_name = get_input(country_i, DATA_V="EXT_data")
        # use AEJ's discount value
        # PFILE = pd.read_csv(f'inputs/AEJ_data/TPFILE', delimiter = "     ", header=None, dtype=float, engine='python').values
        if not country_name == 'NA':
            FUNC_NLLK(PARM=XGUESSES[country_i], CDS1Y=CDS1Y, CDS3Y=CDS3Y, CDS5Y=CDS5Y, PFILE=PFILE, country_name=country_name,
                        XB_FLAG=1, LN_FLAG=1, logging=False, log_file_path="", print_flag=True, CDSLog_flag=False)

    exit(0)

    ###############################################################
    ################# solve parameters via MLE ##################
    ###############################################################
    country_i = 0

    # method = 'L-BFGS-B' # stuck in local min, slow
    # method = 'TNC' # slow
    method = 'SLSQP' # unstable
    # method = 'trust-constr'

    method_options = {'SLSQP':dict(maxiter=5000, disp=True), 
                      'trust-constr': dict(maxiter=5000, verbose=2)}

    now = datetime.now()

    log_file_path = f"logs/{method}/{country_names[country_i]}"
    pathlib.Path(log_file_path).mkdir(parents=True, exist_ok=True)
    log_file_name = f"log_{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}.csv"
    full_log_path = f"{log_file_path}/{log_file_name}"

    with open(full_log_path, "a") as f:
        f.write(f"\tstart time: {now}\n")

    # res = minimize_parallel(FUNC, args=(country_i, True, full_log_path, False, False), x0=XGUESS, bounds=Bounds(lb=XLB, ub=XUB), options=dict(maxiter=5000, ftol=1e-5))
    res = minimize(FUNC_NLLK, args=(country_i, True, full_log_path, False, False), x0=[3.5, -0.8, 0.7, -3.4, 0.6, 0.006, 0.004], method=method, bounds=Bounds(lb=XLB, ub=XUB), options=method_options[method])

    with open(full_log_path, "a") as f:
        f.write(f"\n\tend time: {datetime.now()}\n")

