import os
import pandas as pd




def check_validation(fp):
    sampleid = os.path.basename(fp)
    sampleid = sampleid.replace(".txt", "")
    with open(fp, 'r') as f:
        line = f.readline()
        while line:
            if 'Error found' in line:
                return sampleid, 'Failed'
            elif 'no errors found' in line:
                return sampleid, 'Passed'
            line = f.readline()
        return sampleid, 'No log data found in file'

def main():
    files = snakemake.input
    results = [check_validation(f) for f in files]
    df = pd.DataFrame.from_dict({'patient' : [t[0] for t in results], 
        'result' : [t[1] for t in results]})
    print(df)
    df.to_csv(snakemake.output[0], index = False)



if __name__ == '__main__':
    main()
