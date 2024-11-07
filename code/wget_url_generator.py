# -*- coding: utf-8 -*-
"""
wget scipt url generator
"""

import numpy as np
#%%

r_values=np.int32(np.linspace(1,6,6))
i_values=np.array([1])
p_values=np.array([1])
f_values=np.array([2])

#%%
url_start='http://esgf-data.dkrz.de/esg-search/wget?&download_structure=source_id,experiment_id,variable'
#%%
variable='msftyz'
table_id='Omon'
source_id='CNRM-CM6-1'
    experiment_id='historical'

url= url_start+ '&variable=' + variable+ '&table_id='+table_id+'&source_id='+source_id+'&experiment_id='+experiment_id
#%%Variant labels
for r in range(len(r_values)):
    for i in range (len(i_values)):
        for p in range (len(p_values)):
            for f in range (len(f_values)):
                url=url + '&variant_label=r'+ str(r_values[r]) + 'i' + str(i_values[i])+ 'p' + str(p_values[p])+ 'f' + str(f_values[f])

#%%
    