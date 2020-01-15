import os, shutil

script = os.path.dirname(os.path.realpath("__file__"))
root = os.path.dirname(script)

process_root = root + '/processing'
process_ecmwf = str(process_root + '/process_ecmwf')
process_amsr2 = str(process_root + '/process_amsr2')
def del_file(folder):
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(e)

del_file(process_ecmwf)
del_file(process_amsr2)

