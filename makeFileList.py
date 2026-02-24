import json

failedFilesDict = {}
with open("failedFiles_19Feb.json", "r") as f:
    failed_files = json.load(f)
    failedFilesDict = failed_files
    #print(failed_files)
#print(f"Failed files: {failedFilesDict}")

for key,value in failedFilesDict.items():
    print(f"Key: {key}, Value: {value}")
    with open("fileList_failed19Feb.txt", "a") as f:
        for val in value:
            f.write(f"{val}\n")
