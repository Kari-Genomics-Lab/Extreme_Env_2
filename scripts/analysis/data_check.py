import os
import pandas as pd


DATA_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
RES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results'))


data_temp = pd.read_csv(f"{DATA_PATH}/fragments_100000/Temperature/Extremophiles_Temperature_Summary.tsv", sep='\t')
data_ph = pd.read_csv(f"{DATA_PATH}/fragments_100000/pH/Extremophiles_pH_Summary.tsv", sep='\t')

counter = 0
for i in data_temp.iterrows():
    if i[1]['Assembly'] in list(data_ph['Assembly']):
        # counter+=1
        j = data_ph[data_ph['Assembly'] == i[1]['Assembly']]
        if i[1]['Temperature'] != 'Mesophiles' and i[1]['Domain'] == 'Bacteria':
            counter+=1
            print(i[1]['Assembly'], i[1]['Domain'], i[1]['species'], i[1]['Temperature'], list(j['pH'])[0])

print(counter)
