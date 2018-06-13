import numpy as np
import matplotlib.pyplot as plt

import sys, string, bisect

import pandas as pd

from collections import deque

from scipy.stats import norm


def read_data(filename):

    f = open(filename, 'r')
    data = []

    count = 0

    while True:
        count += 1


        line = f.readline()
        if line == '':
            break
        elements = line.replace('\n','')
        data.append( float(elements) )


    f.close()

    return data


def gen_data(m_size=1000, n_size=128, filename="data.txt"):
    
    index_1 = int(0.2*m_size)
    index_2 = int(0.7*m_size)
    
    motif_test = np.random.normal(50., 10.,size=n_size)
    motif_test_2 = motif_test+np.random.normal(0.,5.,size=len(motif_test))

    data = list(np.random.normal(50., 5.,size=m_size))
    data[index_1:index_1] = motif_test
    data[index_2:index_2] = motif_test_2
    
    f = open(filename, 'w')

    count = 0

    for el in data:
        count += 1
        if count == index_1 + 5:
            el += el
        f.write("%s\n"%el)
        
        
    f.close()



def not_degen_motif(array, epsilon=10.):
    '''
    Check if motif candidate is a degeneritive motif by fitting
    it gainst a threshold.
    '''
    x, y = zip(*array)
    return np.polyfit(x,y,1, full=True)[1] > epsilon
    


def apply_paa(array, w):
    '''
    Piecewise Aggregate Approximation coefs computation.
    '''
    n = len(array)
    c_paa = []

    for i in range(w):
        nw = int(1.0*n/w)
        intern_sum_range = range(nw*i,(i+1)*nw)


        suma = 0.0
        for j in intern_sum_range:
            suma += array[j]
        
        c_paa.append(suma/nw)

    return c_paa


def gen_breakpoints(a):
    result = []
    cdf_beta_old = 0
    for index in range(1,a):

        beta_new = norm.ppf(1./a + cdf_beta_old)
        cdf_beta_old = norm.cdf(beta_new)
        result.append(beta_new)
    
    return result


def apply_symb(coefs, symbs):

    a = len(symbs)
    result = []
    norm_coefs = normalize_cs(coefs)
    betas = gen_breakpoints(a)

    for c in norm_coefs:
        index = bisect.bisect_left(betas, c)
        result.append(symbs[index])

    return result



def normalize_cs(coefs):
    sd = np.std(coefs)
    mean = np.mean(coefs)
    return (np.array(coefs) - mean)/sd



def df_get_col_rle(df, i):

    if i in df.columns.values.tolist():
        return df[i]
    else:
        return df_get_col_rle(df, i-1)



class TSP1(object):


    def __init__(self, data, n, w, a):
        
        self.m = len(data)
        self.data = np.array(data)
        self.n = n
        self.w = w
        self.symbs = list(string.ascii_lowercase)[:a]

        self.paa_matrix , self.s_matrix= self.treat_data()
        


    def out_s1(self, ind):
        
        self.plot_window(ind, self.paa_matrix[ind], filename='test_window_%s.png'%ind, \
            tit_str=''.join(df_get_col_rle(self.s_matrix, ind).tolist()))

        self.s_matrix.T.to_csv('S_hat.txt', sep='\t', encoding='utf-8')


    def treat_data(self):

        c_bar = {}
        c_hat = {}

        prev = [0.]*self.w
        comp_count = 0

        for i in range(self.m-self.n+1):
            window_i = (self.data[i:i+self.n])
            
            # APPLY LINEARITY FILTER
            

            c_bar[i] = apply_paa(window_i, self.w)

            curr = apply_symb(c_bar[i], self.symbs)        

            
            # RLE compression
            if prev == curr:
                comp_count += 1
            else:
                c_hat[i] = curr
                prev = c_hat[i]
        

            #c_hat[i] = curr (for no compression)
                
        paa_matrix = pd.DataFrame(c_bar, index=range(self.w))

        s_matrix = pd.DataFrame(c_hat, index=range(self.w))


        print "RLE Compression: %.1f%%" %((100.*comp_count)/(self.m-self.n+1))
        print "-----------------------------------------------"

        return paa_matrix, s_matrix




    def plot_window(self, index, c_bar, filename='test_window.png', tit_str=''):

        fig = plt.figure(figsize=(10,10))

        ax = fig.add_subplot(111)

        ind = range(index,index+self.n)
        d = self.data[index:index+self.n]


        mean = np.mean(c_bar)
        sd = np.std(c_bar)

        d = (np.array(d)-mean)/sd

        plt.plot(ind, d, color='#1e90ff', marker='', linestyle='-', lw=2.5, mew=1.)

        count = 0
        
        c_bar = (np.array(c_bar)-mean)/sd
        
        for el in c_bar:
            c_bar_plot = []
            for i in range(int(self.n/self.w)):
                c_bar_plot.append( (index+count, el) )
                count += 1

            plt.plot(*zip(*c_bar_plot), color='#ff4500', marker='', linestyle='-', lw=3.5, mew=1.)


        betas = gen_breakpoints(len(self.symbs))

        for b in betas:
            plt.axhline(b, color='k', ls='--', lw = 2.)
            

        plt.xlabel('t')
        plt.ylabel('x')
        plt.suptitle(tit_str+'; index=%s, n=%s; w=%s'%(index, self.n, self.w))
        fig.savefig(filename, format='png')
        plt.close()


if __name__ == '__main__':

    gen_data()
    

































