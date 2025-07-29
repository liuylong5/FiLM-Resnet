import os
import shutil
import time
import glob


BaseNNDir = ""
BaseBWDDir = ""
WorkDir = ""
TmpDir = ""

def initialize_evaluation(RawModel):
    global BaseNNDir, BaseBWDDir, BaseFWDDir, WorkDir, TmpDir 
    WorkDir = os.getcwd()
    UtilDir = '../jacobian/standare_batch_autograd/util_new'
    RawBWD = '../jacobian/standare_batch_bwd'
    RawFWD = '../jacobian/standare_batch_fwd'

    BaseNNDir = os.path.join(WorkDir, 'nn')
    BaseBWDDir = os.path.join(WorkDir, 'bwd')
    BaseFWDDir = os.path.join(WorkDir, 'fwd')

    src1 = os.path.join(RawModel, 'src')
    dst1 = os.path.join(BaseNNDir, 'src')
    src2 = UtilDir
    dst2 = os.path.join(BaseNNDir, 'util')
    src3 = os.path.join(RawBWD, 'src_raw')
    dst3 = os.path.join(BaseBWDDir, 'src_raw')
    dst4 = os.path.join(BaseBWDDir, 'src')

    src4 = os.path.join(RawFWD, 'src_raw')
    dst5 = os.path.join(BaseFWDDir, 'src_raw')
    dst6 = os.path.join(BaseFWDDir, 'src')

    folders = [
        (src1, dst1),
        (src2, dst2),
        (src3, dst3),
        (src3, dst4),
        (src4, dst5),
        (src4, dst6)
    ]

    for src, dst in folders:
        if os.path.exists(dst):
            shutil.rmtree(dst) 
        shutil.copytree(src, dst)  
    shutil.copy(os.path.join(RawBWD, 'makefile.intel'), os.path.join(BaseBWDDir, 'makefile.intel'))
    shutil.copy(os.path.join(RawBWD, 'path.csh'), os.path.join(BaseBWDDir, 'path.csh'))
    shutil.copy(os.path.join(RawFWD, 'makefile.intel'), os.path.join(BaseFWDDir, 'makefile.intel'))
    shutil.copy(os.path.join(RawFWD, 'path.csh'), os.path.join(BaseFWDDir, 'path.csh'))
    TmpDir = os.path.join(WorkDir, 'tmp')
    os.makedirs(TmpDir, exist_ok=True)
    print(f"{RawModel = }")
    print(f"{BaseNNDir = }")
    print(f"{BaseBWDDir = }")
    print(f"{BaseFWDDir = }")
    print(f"{WorkDir = }")
    return 'initialize evaluation done'

def wait_for_log(log_file, target_string, timeout=600, interval=5):
    start_time = time.time()
    while time.time() - start_time < timeout:
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                if target_string in f.read():
                    return True
        time.sleep(interval)
    raise TimeoutError(f"in: {log_file} no '{target_string}'")

def get_raw_data(h5_file, save_size=10):
    data_to_npy_str = ""
    data_to_npy_str += "import sys\n"
    data_to_npy_str += f"sys.path.insert(0, '{BaseNNDir}')\n"  
    with open(f"{BaseNNDir}/util/data_to_npy.py", 'r') as f:
        for line in f:
            if "save_size =" in line:
                data_to_npy_str += f"save_size = {save_size}\n"
            elif "hdf5_file_lists =" in line:
                data_to_npy_str += f"hdf5_file_lists = ['{h5_file}']\n"
            elif "params_array =" in line:
                data_to_npy_str += line
                data_to_npy_str += f"conc_72 = numpy.full((inputs_array.shape[0], 72), 1e-19)\n"
                data_to_npy_str += f"remaining_idx_lists = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]\n"
                data_to_npy_str += f"conc_72[:, remaining_idx_lists] = inputs_array\n"
                data_to_npy_str += f"inputs_array = conc_72\n"
            elif "numpy.save('inputs.npy', inputs_array)" in line:
                data_to_npy_str += f"numpy.save(os.path.join('{TmpDir}', 'inputs.npy'), inputs_array)\n"
            elif "numpy.save('params.npy', params_array)" in line:
                data_to_npy_str += f"numpy.save(os.path.join('{TmpDir}', 'params.npy'), params_array)\n"
            else:
                data_to_npy_str += line
    with open(f"{TmpDir}/data_to_npy.py", 'w') as f:
        f.write(data_to_npy_str)
    data_to_npy_log = os.path.join(f'{TmpDir}', 'data_to_npy.log')
    if os.path.exists(data_to_npy_log): 
        os.remove(data_to_npy_log) 
    os.system(f"/dssg/home/acct-esehazenet/hazenet-liuyulong/.conda/envs/pytorch-env/bin/python -u {TmpDir}/data_to_npy.py >> {data_to_npy_log}")
    wait_for_log(data_to_npy_log, 'Process ends at')
    return 'get_raw_data done'

def data_to_csv(): 
    data_to_csv_str = ""
    with open(f"{BaseNNDir}/util/data_to_csv.py", 'r') as f:
        for line in f:
            if "inputs_npy =" in line:
                data_to_csv_str += f"inputs_npy = os.path.join('{TmpDir}', 'inputs.npy')\n"
            elif "inputs_csv =" in line:
                data_to_csv_str += f"inputs_csv = os.path.join('{TmpDir}', 'inputs.csv')\n"
            elif "params_npy =" in line:
                data_to_csv_str += f"params_npy = os.path.join('{TmpDir}', 'params.npy')\n"
            elif "params_csv =" in line:
                data_to_csv_str += f"params_csv = os.path.join('{TmpDir}', 'params.csv')\n"
            else:
                data_to_csv_str += line
    with open(f"{TmpDir}/data_to_csv.py", 'w') as f:
        f.write(data_to_csv_str)
    data_to_csv_log = os.path.join(f'{TmpDir}', 'data_to_csv.log')
    if os.path.exists(data_to_csv_log): 
        os.remove(data_to_csv_log) 
    os.system(f"/dssg/home/acct-esehazenet/hazenet-liuyulong/.conda/envs/pytorch-env/bin/python -u {TmpDir}/data_to_csv.py >> {data_to_csv_log}")
    wait_for_log(data_to_csv_log, 'CSV array saved')
    return 'data_to_csv done'

def modify_bwd(save_size=10):
    modify_bwd_str = ""
    with open(f"{BaseBWDDir}/src_raw/kppinit.F", 'r') as f:
        for line in f:
            if "OPEN(UNIT=10086" in line:
                file_name = os.path.join(TmpDir, 'inputs.csv')
                modify_bwd_str += f"      BASEFILE = trim('{file_name[:44]}') // \n"
                modify_bwd_str += f"     &  trim('{file_name[44:]}')\n\n"
                modify_bwd_str += f'      OPEN(UNIT=10086, FILE=BASEFILE,\n'
                modify_bwd_str += f'     &  STATUS="OLD", ACTION="READ", IOSTAT=iostat)\n'
            elif "OPEN(UNIT=10087" in line:
                file_name = os.path.join(TmpDir, 'params.csv')
                modify_bwd_str += f"      BASEFILE = trim('{file_name[:44]}') // \n"
                modify_bwd_str += f"     &  trim('{file_name[44:]}')\n\n"
                modify_bwd_str += f'      OPEN(UNIT=10087, FILE=BASEFILE,\n'
                modify_bwd_str += f'     &  STATUS="OLD", ACTION="READ", IOSTAT=iostat)\n'
            else:
                modify_bwd_str += line
    with open(f"{BaseBWDDir}/src/kppinit.F", 'w') as f:
        f.write(modify_bwd_str)

    modify_bwd_str = ""
    with open(f"{BaseBWDDir}/src_raw/kppdata_mod.F", 'r') as f:
        for line in f:
            if "NUM_CASES =" in line:
                modify_bwd_str += f'      INTEGER, PARAMETER :: NUM_CASES = {save_size}\n'
            else:
                modify_bwd_str += line
    with open(f"{BaseBWDDir}/src/kppdata_mod.F", 'w') as f:
        f.write(modify_bwd_str)

    compile_bwd_log = os.path.join(f'{TmpDir}', 'compile_bwd.log')
    run_bwd_log = os.path.join(f'{TmpDir}', 'run_bwd.log')
    for log_file in [compile_bwd_log, run_bwd_log]:
        if os.path.exists(log_file): 
            os.remove(log_file) 

    info_path = '../jacobian/standare_batch_bwd/path.csh'

    os.chdir(BaseBWDDir)

    file_patterns = ["*.o", "*.mod", "bwd_intel"]
    files_to_delete = []
    for pattern in file_patterns:
        files_to_delete.extend(glob.glob(pattern))  

    for file in files_to_delete:
        file_path = os.path.join(BaseBWDDir, file)
        if os.path.exists(file_path):
            os.remove(file_path)
        else:
            print(f"no files: {file_path}")

    os.system(f"cd {BaseBWDDir} && module purge && module load oneapi && source {info_path} && make -f makefile.intel bwd_intel >> '{compile_bwd_log}'")

    wait_for_log(compile_bwd_log, '-lmpifort -qopenmp -o bwd_intel')

    os.system(f"cd {BaseBWDDir} && module purge && module load oneapi && source {info_path}  && ./bwd_intel >> '{run_bwd_log}'")

    wait_for_log(run_bwd_log, 'jacobians.csv')

    source_jacobian = os.path.join(BaseBWDDir, 'jacobians.csv')
    destination_jacobian = os.path.join(TmpDir, 'jacobian.csv')
    shutil.copy(source_jacobian, destination_jacobian)
    return 'modify and recompile bwd done'


def save_bwd_jacobian(save_size=10):
    save_bwd_jacobian_str = ""
    with open(f"{BaseNNDir}/util/save_bwd_jacobian.py", 'r') as f:
        for line in f:
            if "save_size =" in line:
                save_bwd_jacobian_str += f"save_size = {save_size}\n"
            elif "save_path =" in line:

                save_bwd_jacobian_str += f"save_path = os.path.join('{WorkDir}', 'bwd_jacobian.npy')\n"
            elif "np.save(save_path, jacobians)" in line:
                save_bwd_jacobian_str += line
                save_bwd_jacobian_str += f"remaining_idx_lists = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 54, 55, 66, 67, 70]\n"
                save_bwd_jacobian_str += f"J_new = jacobians[:, remaining_idx_lists, :][:, :, remaining_idx_lists]\n"
                save_bwd_jacobian_str += f"save_path = os.path.join('{WorkDir}', 'bwd_jacobian.npy')\n"
                save_bwd_jacobian_str += f"np.save(save_path, J_new)\n"
            else:
                save_bwd_jacobian_str += line
    with open(f"{TmpDir}/save_bwd_jacobian.py", 'w') as f:
        f.write(save_bwd_jacobian_str)
    save_bwd_jacobian_log = os.path.join(f'{TmpDir}', 'save_bwd_jacobian.log')
    if os.path.exists(save_bwd_jacobian_log): 
        os.remove(save_bwd_jacobian_log) 
    os.system(f"/dssg/home/acct-esehazenet/hazenet-liuyulong/.conda/envs/pytorch-env/bin/python -u {TmpDir}/save_bwd_jacobian.py >> {save_bwd_jacobian_log}")
    wait_for_log(save_bwd_jacobian_log, 'Jacobian saved!')
    return 'save bwd npy done'

def modify_fwd(save_size=10):
    modify_fwd_str = ""
    with open(f"{BaseFWDDir}/src_raw/kppinit.F", 'r') as f:
        for line in f:
            if "OPEN(UNIT=10086" in line:
                file_name = os.path.join(TmpDir, 'inputs.csv')
                modify_fwd_str += f"      BASEFILE = trim('{file_name[:44]}') // \n"
                modify_fwd_str += f"     &  trim('{file_name[44:]}')\n\n"
                modify_fwd_str += f'      OPEN(UNIT=10086, FILE=BASEFILE,\n'
                modify_fwd_str += f'     &  STATUS="OLD", ACTION="READ", IOSTAT=iostat)\n'
            elif "OPEN(UNIT=10087" in line:
                file_name = os.path.join(TmpDir, 'params.csv')
                modify_fwd_str += f"      BASEFILE = trim('{file_name[:44]}') // \n"
                modify_fwd_str += f"     &  trim('{file_name[44:]}')\n\n"
                modify_fwd_str += f'      OPEN(UNIT=10087, FILE=BASEFILE,\n'
                modify_fwd_str += f'     &  STATUS="OLD", ACTION="READ", IOSTAT=iostat)\n'
            else:
                modify_fwd_str += line
    with open(f"{BaseFWDDir}/src/kppinit.F", 'w') as f:
        f.write(modify_fwd_str)

    modify_fwd_str = ""
    with open(f"{BaseFWDDir}/src_raw/kppdata_mod.F", 'r') as f:
        for line in f:
            if "NUM_CASES =" in line:
                modify_fwd_str += f'      INTEGER, PARAMETER :: NUM_CASES = {save_size}\n'
            else:
                modify_fwd_str += line
    with open(f"{BaseFWDDir}/src/kppdata_mod.F", 'w') as f:
        f.write(modify_fwd_str)

    compile_fwd_log = os.path.join(f'{TmpDir}', 'compile_fwd.log')
    run_fwd_log = os.path.join(f'{TmpDir}', 'run_fwd.log')
    for log_file in [compile_fwd_log, run_fwd_log]:
        if os.path.exists(log_file): 
            os.remove(log_file) 

    info_path = '../jacobian/standare_batch_fwd/path.csh'

    os.chdir(BaseFWDDir)

    file_patterns = ["*.o", "*.mod", "main"]
    files_to_delete = []
    for pattern in file_patterns:
        files_to_delete.extend(glob.glob(pattern))  

    for file in files_to_delete:
        file_path = os.path.join(BaseFWDDir, file)
        if os.path.exists(file_path):
            os.remove(file_path)
        else:
            print(f"no files: {file_path}")

    os.system(f"cd {BaseFWDDir} && module purge && module load oneapi && source {info_path} && make -f makefile.intel main >> '{compile_fwd_log}'")

    wait_for_log(compile_fwd_log, '-lmpi -qopenmp -o main')

    os.system(f"cd {BaseFWDDir} && module purge && module load oneapi && source {info_path}  && ./main >> '{run_fwd_log}'")

    wait_for_log(run_fwd_log, 'jacobians.csv')

    source_jacobian = os.path.join(BaseFWDDir, 'jacobians.csv')
    destination_jacobian = os.path.join(TmpDir, 'jacobian.csv')
    shutil.copy(source_jacobian, destination_jacobian)
    return 'modify and recompile fwd done'

def save_fwd_jacobian(save_size=10):
    save_fwd_jacobian_str = ""
    with open(f"{BaseNNDir}/util/save_fwd_jacobian.py", 'r') as f:
        for line in f:
            if "save_size =" in line:
                save_fwd_jacobian_str += f"save_size = {save_size}\n"
            elif "save_path =" in line:
                save_fwd_jacobian_str += f"save_path = os.path.join('{WorkDir}', 'fwd_jacobian.npy')\n"
            elif "np.save(save_path, jacobians)" in line:
                save_fwd_jacobian_str += line
                save_fwd_jacobian_str += f"remaining_idx_lists = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 49, 51, 52, 54, 55, 66, 67, 70]\n"
                save_fwd_jacobian_str += f"J_new = jacobians[:, remaining_idx_lists, :][:, :, remaining_idx_lists]\n"
                save_fwd_jacobian_str += f"save_path = os.path.join('{WorkDir}', 'fwd_jacobian.npy')\n"
                save_fwd_jacobian_str += f"np.save(save_path, J_new)\n"
            else:
                save_fwd_jacobian_str += line
    with open(f"{TmpDir}/save_fwd_jacobian.py", 'w') as f:
        f.write(save_fwd_jacobian_str)
    save_fwd_jacobian_log = os.path.join(f'{TmpDir}', 'save_fwd_jacobian.log')
    if os.path.exists(save_fwd_jacobian_log): 
        os.remove(save_fwd_jacobian_log) 
    os.system(f"/dssg/home/acct-esehazenet/hazenet-liuyulong/.conda/envs/pytorch-env/bin/python -u {TmpDir}/save_fwd_jacobian.py >> {save_fwd_jacobian_log}")
    wait_for_log(save_fwd_jacobian_log, 'Jacobian saved!')
    return 'save fwd npy done'


def save_nn_jacobian(model_path):
    save_nn_jacobian_str = ""
    save_nn_jacobian_str += "import sys\n"
    save_nn_jacobian_str += f"sys.path.insert(0, '{BaseNNDir}')\n" 
    with open(f"{BaseNNDir}/util/save_nn_jacobian.py", 'r') as f:
        for line in f:
            if "inputs_npy =" in line:
                save_nn_jacobian_str += f"inputs_npy = os.path.join('{TmpDir}', 'inputs.npy')\n"
            elif "params_npy =" in line:
                save_nn_jacobian_str += f"params_npy = os.path.join('{TmpDir}', 'params.npy')\n"
            elif "model_path =" in line:
                save_nn_jacobian_str += f"model_path = '{model_path}'\n"
            elif "norm_path =" in line:
                save_nn_jacobian_str += f"norm_path = os.path.join('{model_path.split('path')[0]}', 'data/norm_values.pkl')\n"
            elif "nn_jacobian_path =" in line:
                save_nn_jacobian_str += f"nn_jacobian_path = os.path.join('{WorkDir}', 'nn_jacobian.npy')\n"
            else:
                save_nn_jacobian_str += line
    with open(f"{TmpDir}/save_nn_jacobian.py", 'w') as f:
        f.write(save_nn_jacobian_str)
    save_nn_jacobian_log = os.path.join(f'{TmpDir}', 'save_nn_jacobian.log')
    if os.path.exists(save_nn_jacobian_log): 
        os.remove(save_nn_jacobian_log) 
    os.system(f"/dssg/home/acct-esehazenet/hazenet-liuyulong/.conda/envs/pytorch-env/bin/python -u {TmpDir}/save_nn_jacobian.py >> {save_nn_jacobian_log}")
    wait_for_log(save_nn_jacobian_log, 'Process ends at')
    return 'save_nn_jacobian done'

def main():
    RawModel = 'model'
    model_name = '20250617-231942.pt'
    h5_file = '../merged_datasets/test_dataset_001.h5'
    model_path = os.path.join(RawModel, 'path', model_name)
    save_size = 1000

    tasks = [
        (initialize_evaluation, [RawModel]),
        (get_raw_data, [h5_file, save_size]),
        (data_to_csv, []),
        (modify_fwd, [save_size]),
        (save_fwd_jacobian, [save_size]),      
    ]

    for func, args in tasks:
        message = func(*args)
        print(message)

    print('all done!')


if __name__ == "__main__":
    main()