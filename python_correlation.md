```python
#Import the packages
from pandas import *
import numpy as np
from scipy.stats.stats import pearsonr
import itertools
```


```python
## Create a dataframe as like as gene expression matrix
df = DataFrame(np.random.random((10, 10)), columns=['gene_' + chr(i + ord('a')) for i in range(10)],
              index=['Sample_' + chr(i + ord('k')) for i in range(10)]) 
df
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_a</th>
      <th>gene_b</th>
      <th>gene_c</th>
      <th>gene_d</th>
      <th>gene_e</th>
      <th>gene_f</th>
      <th>gene_g</th>
      <th>gene_h</th>
      <th>gene_i</th>
      <th>gene_j</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Sample_k</th>
      <td>0.318694</td>
      <td>0.876427</td>
      <td>0.207955</td>
      <td>0.900428</td>
      <td>0.903615</td>
      <td>0.564069</td>
      <td>0.196870</td>
      <td>0.695175</td>
      <td>0.259064</td>
      <td>0.839468</td>
    </tr>
    <tr>
      <th>Sample_l</th>
      <td>0.686489</td>
      <td>0.910583</td>
      <td>0.698854</td>
      <td>0.941880</td>
      <td>0.141722</td>
      <td>0.299639</td>
      <td>0.645385</td>
      <td>0.843439</td>
      <td>0.153848</td>
      <td>0.470428</td>
    </tr>
    <tr>
      <th>Sample_m</th>
      <td>0.490361</td>
      <td>0.916693</td>
      <td>0.921328</td>
      <td>0.981796</td>
      <td>0.300647</td>
      <td>0.191906</td>
      <td>0.362133</td>
      <td>0.049860</td>
      <td>0.599518</td>
      <td>0.962537</td>
    </tr>
    <tr>
      <th>Sample_n</th>
      <td>0.219756</td>
      <td>0.187282</td>
      <td>0.556514</td>
      <td>0.107674</td>
      <td>0.480256</td>
      <td>0.493407</td>
      <td>0.316614</td>
      <td>0.528466</td>
      <td>0.067806</td>
      <td>0.357675</td>
    </tr>
    <tr>
      <th>Sample_o</th>
      <td>0.639952</td>
      <td>0.148242</td>
      <td>0.288024</td>
      <td>0.888205</td>
      <td>0.538402</td>
      <td>0.983597</td>
      <td>0.179698</td>
      <td>0.299926</td>
      <td>0.899446</td>
      <td>0.375380</td>
    </tr>
    <tr>
      <th>Sample_p</th>
      <td>0.035375</td>
      <td>0.848469</td>
      <td>0.290245</td>
      <td>0.628137</td>
      <td>0.476790</td>
      <td>0.566463</td>
      <td>0.183813</td>
      <td>0.774419</td>
      <td>0.393618</td>
      <td>0.670081</td>
    </tr>
    <tr>
      <th>Sample_q</th>
      <td>0.984982</td>
      <td>0.453418</td>
      <td>0.714330</td>
      <td>0.673423</td>
      <td>0.517886</td>
      <td>0.951160</td>
      <td>0.363136</td>
      <td>0.360941</td>
      <td>0.599349</td>
      <td>0.162177</td>
    </tr>
    <tr>
      <th>Sample_r</th>
      <td>0.747231</td>
      <td>0.409106</td>
      <td>0.221303</td>
      <td>0.864283</td>
      <td>0.518088</td>
      <td>0.712261</td>
      <td>0.689450</td>
      <td>0.038962</td>
      <td>0.623002</td>
      <td>0.176004</td>
    </tr>
    <tr>
      <th>Sample_s</th>
      <td>0.530720</td>
      <td>0.491937</td>
      <td>0.634272</td>
      <td>0.337528</td>
      <td>0.375218</td>
      <td>0.356336</td>
      <td>0.498910</td>
      <td>0.924223</td>
      <td>0.584844</td>
      <td>0.243460</td>
    </tr>
    <tr>
      <th>Sample_t</th>
      <td>0.882852</td>
      <td>0.748447</td>
      <td>0.601405</td>
      <td>0.778175</td>
      <td>0.616356</td>
      <td>0.247548</td>
      <td>0.180486</td>
      <td>0.625956</td>
      <td>0.611697</td>
      <td>0.703179</td>
    </tr>
  </tbody>
</table>
</div>




```python
correlations = {}
#store the gene names
columns = df.columns.tolist()

#perform pearson correlation
for col_a, col_b in itertools.combinations(columns, 2):
    correlations[col_a +':' + col_b] = pearsonr(df.loc[:, col_a], df.loc[:, col_b])

#Convert the dictionary into dataframe
result = DataFrame.from_dict(correlations, orient='index')

#Assign the column names
result.columns = ['PCC', 'p-value']

print(result.sort_index())
```

                        PCC   p-value
    gene_a:gene_b -0.153370  0.672284
    gene_a:gene_c  0.316841  0.372398
    gene_a:gene_d  0.356856  0.311421
    gene_a:gene_e -0.103400  0.776216
    gene_a:gene_f  0.192147  0.594857
    gene_a:gene_g  0.328188  0.354544
    gene_a:gene_h -0.304845  0.391742
    gene_a:gene_i  0.497258  0.143658
    gene_a:gene_j -0.395864  0.257467
    gene_b:gene_c  0.238075  0.507739
    gene_b:gene_d  0.501916  0.139352
    gene_b:gene_e -0.098959  0.785634
    gene_b:gene_f -0.628455  0.051650
    gene_b:gene_g -0.010050  0.978018
    gene_b:gene_h  0.289964  0.416399
    gene_b:gene_i -0.294797  0.408312
    gene_b:gene_j  0.733917  0.015675
    gene_c:gene_d -0.068959  0.849868
    gene_c:gene_e -0.635387  0.048363
    gene_c:gene_f -0.515245  0.127472
    gene_c:gene_g  0.204358  0.571174
    gene_c:gene_h -0.036244  0.920821
    gene_c:gene_i -0.054534  0.881060
    gene_c:gene_j  0.123117  0.734726
    gene_d:gene_e  0.032626  0.928707
    gene_d:gene_f  0.026852  0.941305
    gene_d:gene_g  0.062359  0.864118
    gene_d:gene_h -0.364520  0.300385
    gene_d:gene_i  0.354651  0.314635
    gene_d:gene_j  0.419242  0.227813
    gene_e:gene_f  0.361560  0.304622
    gene_e:gene_g -0.583162  0.076806
    gene_e:gene_h -0.028750  0.937161
    gene_e:gene_i  0.080331  0.825406
    gene_e:gene_j  0.181629  0.615535
    gene_f:gene_g -0.146843  0.685619
    gene_f:gene_h -0.321660  0.364761
    gene_f:gene_i  0.391805  0.262821
    gene_f:gene_j -0.560116  0.092192
    gene_g:gene_h -0.117758  0.745946
    gene_g:gene_i -0.137741  0.704343
    gene_g:gene_j -0.497628  0.143313
    gene_h:gene_i -0.481393  0.158926
    gene_h:gene_j  0.036145  0.921035
    gene_i:gene_j -0.161884  0.655014
    


```python
#Convert to a dataframe
result = DataFrame.from_dict(correlations, orient='index')
result.columns = ['PCC', 'p-value']
#Convert index as a column
result.reset_index(level=0, inplace=True)
#rename the converted column
result.rename(columns = {'index':'Genes'}, inplace = True)
#split the column as two columns
result[['Gene-1','Gene-2']] = result.Genes.str.split(":", expand=True)
#select the required variables
result = result[['Gene-1', 'Gene-2', 'PCC', 'p-value']]
result
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene-1</th>
      <th>Gene-2</th>
      <th>PCC</th>
      <th>p-value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>gene_a</td>
      <td>gene_b</td>
      <td>-0.153370</td>
      <td>0.672284</td>
    </tr>
    <tr>
      <th>1</th>
      <td>gene_a</td>
      <td>gene_c</td>
      <td>0.316841</td>
      <td>0.372398</td>
    </tr>
    <tr>
      <th>2</th>
      <td>gene_a</td>
      <td>gene_d</td>
      <td>0.356856</td>
      <td>0.311421</td>
    </tr>
    <tr>
      <th>3</th>
      <td>gene_a</td>
      <td>gene_e</td>
      <td>-0.103400</td>
      <td>0.776216</td>
    </tr>
    <tr>
      <th>4</th>
      <td>gene_a</td>
      <td>gene_f</td>
      <td>0.192147</td>
      <td>0.594857</td>
    </tr>
    <tr>
      <th>5</th>
      <td>gene_a</td>
      <td>gene_g</td>
      <td>0.328188</td>
      <td>0.354544</td>
    </tr>
    <tr>
      <th>6</th>
      <td>gene_a</td>
      <td>gene_h</td>
      <td>-0.304845</td>
      <td>0.391742</td>
    </tr>
    <tr>
      <th>7</th>
      <td>gene_a</td>
      <td>gene_i</td>
      <td>0.497258</td>
      <td>0.143658</td>
    </tr>
    <tr>
      <th>8</th>
      <td>gene_a</td>
      <td>gene_j</td>
      <td>-0.395864</td>
      <td>0.257467</td>
    </tr>
    <tr>
      <th>9</th>
      <td>gene_b</td>
      <td>gene_c</td>
      <td>0.238075</td>
      <td>0.507739</td>
    </tr>
    <tr>
      <th>10</th>
      <td>gene_b</td>
      <td>gene_d</td>
      <td>0.501916</td>
      <td>0.139352</td>
    </tr>
    <tr>
      <th>11</th>
      <td>gene_b</td>
      <td>gene_e</td>
      <td>-0.098959</td>
      <td>0.785634</td>
    </tr>
    <tr>
      <th>12</th>
      <td>gene_b</td>
      <td>gene_f</td>
      <td>-0.628455</td>
      <td>0.051650</td>
    </tr>
    <tr>
      <th>13</th>
      <td>gene_b</td>
      <td>gene_g</td>
      <td>-0.010050</td>
      <td>0.978018</td>
    </tr>
    <tr>
      <th>14</th>
      <td>gene_b</td>
      <td>gene_h</td>
      <td>0.289964</td>
      <td>0.416399</td>
    </tr>
    <tr>
      <th>15</th>
      <td>gene_b</td>
      <td>gene_i</td>
      <td>-0.294797</td>
      <td>0.408312</td>
    </tr>
    <tr>
      <th>16</th>
      <td>gene_b</td>
      <td>gene_j</td>
      <td>0.733917</td>
      <td>0.015675</td>
    </tr>
    <tr>
      <th>17</th>
      <td>gene_c</td>
      <td>gene_d</td>
      <td>-0.068959</td>
      <td>0.849868</td>
    </tr>
    <tr>
      <th>18</th>
      <td>gene_c</td>
      <td>gene_e</td>
      <td>-0.635387</td>
      <td>0.048363</td>
    </tr>
    <tr>
      <th>19</th>
      <td>gene_c</td>
      <td>gene_f</td>
      <td>-0.515245</td>
      <td>0.127472</td>
    </tr>
    <tr>
      <th>20</th>
      <td>gene_c</td>
      <td>gene_g</td>
      <td>0.204358</td>
      <td>0.571174</td>
    </tr>
    <tr>
      <th>21</th>
      <td>gene_c</td>
      <td>gene_h</td>
      <td>-0.036244</td>
      <td>0.920821</td>
    </tr>
    <tr>
      <th>22</th>
      <td>gene_c</td>
      <td>gene_i</td>
      <td>-0.054534</td>
      <td>0.881060</td>
    </tr>
    <tr>
      <th>23</th>
      <td>gene_c</td>
      <td>gene_j</td>
      <td>0.123117</td>
      <td>0.734726</td>
    </tr>
    <tr>
      <th>24</th>
      <td>gene_d</td>
      <td>gene_e</td>
      <td>0.032626</td>
      <td>0.928707</td>
    </tr>
    <tr>
      <th>25</th>
      <td>gene_d</td>
      <td>gene_f</td>
      <td>0.026852</td>
      <td>0.941305</td>
    </tr>
    <tr>
      <th>26</th>
      <td>gene_d</td>
      <td>gene_g</td>
      <td>0.062359</td>
      <td>0.864118</td>
    </tr>
    <tr>
      <th>27</th>
      <td>gene_d</td>
      <td>gene_h</td>
      <td>-0.364520</td>
      <td>0.300385</td>
    </tr>
    <tr>
      <th>28</th>
      <td>gene_d</td>
      <td>gene_i</td>
      <td>0.354651</td>
      <td>0.314635</td>
    </tr>
    <tr>
      <th>29</th>
      <td>gene_d</td>
      <td>gene_j</td>
      <td>0.419242</td>
      <td>0.227813</td>
    </tr>
    <tr>
      <th>30</th>
      <td>gene_e</td>
      <td>gene_f</td>
      <td>0.361560</td>
      <td>0.304622</td>
    </tr>
    <tr>
      <th>31</th>
      <td>gene_e</td>
      <td>gene_g</td>
      <td>-0.583162</td>
      <td>0.076806</td>
    </tr>
    <tr>
      <th>32</th>
      <td>gene_e</td>
      <td>gene_h</td>
      <td>-0.028750</td>
      <td>0.937161</td>
    </tr>
    <tr>
      <th>33</th>
      <td>gene_e</td>
      <td>gene_i</td>
      <td>0.080331</td>
      <td>0.825406</td>
    </tr>
    <tr>
      <th>34</th>
      <td>gene_e</td>
      <td>gene_j</td>
      <td>0.181629</td>
      <td>0.615535</td>
    </tr>
    <tr>
      <th>35</th>
      <td>gene_f</td>
      <td>gene_g</td>
      <td>-0.146843</td>
      <td>0.685619</td>
    </tr>
    <tr>
      <th>36</th>
      <td>gene_f</td>
      <td>gene_h</td>
      <td>-0.321660</td>
      <td>0.364761</td>
    </tr>
    <tr>
      <th>37</th>
      <td>gene_f</td>
      <td>gene_i</td>
      <td>0.391805</td>
      <td>0.262821</td>
    </tr>
    <tr>
      <th>38</th>
      <td>gene_f</td>
      <td>gene_j</td>
      <td>-0.560116</td>
      <td>0.092192</td>
    </tr>
    <tr>
      <th>39</th>
      <td>gene_g</td>
      <td>gene_h</td>
      <td>-0.117758</td>
      <td>0.745946</td>
    </tr>
    <tr>
      <th>40</th>
      <td>gene_g</td>
      <td>gene_i</td>
      <td>-0.137741</td>
      <td>0.704343</td>
    </tr>
    <tr>
      <th>41</th>
      <td>gene_g</td>
      <td>gene_j</td>
      <td>-0.497628</td>
      <td>0.143313</td>
    </tr>
    <tr>
      <th>42</th>
      <td>gene_h</td>
      <td>gene_i</td>
      <td>-0.481393</td>
      <td>0.158926</td>
    </tr>
    <tr>
      <th>43</th>
      <td>gene_h</td>
      <td>gene_j</td>
      <td>0.036145</td>
      <td>0.921035</td>
    </tr>
    <tr>
      <th>44</th>
      <td>gene_i</td>
      <td>gene_j</td>
      <td>-0.161884</td>
      <td>0.655014</td>
    </tr>
  </tbody>
</table>
</div>




```python

```


```python

```


```python

```
