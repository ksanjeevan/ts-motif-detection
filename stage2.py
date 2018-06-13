from stage1 import *
import random, math





def nCr(n, r):
    '''
    Binomial coefficient computation.
    '''
    f = math.factorial
    return f(n) / (f(r) * f(n-r) )


def expect_s_thresh(num_words, a, w, d, mask_size):
    '''
    Threshold expectation for a single collision matrix entry per iteration, E(k,a,w,d,t).
    (buggy)
    '''
    return sum([((1.-float(i)/w)**mask_size) * nCr(w, i) * ((1.-1./a)**i) \
        * ((1./a)**(w-i)) for i in range(d+1)])*nCr(num_words, 2)


def dist_function(c, m, d):
    a = normalize_cs(np.array(c))
    b = normalize_cs(np.array(m))
    u = sorted(np.absolute(a-b), reverse=True)[d:]

    return np.linalg.norm(u)



class TSP2(TSP1):


    def __init__(self, data, params, k, d):
        
        n, w, a = params
        TSP1.__init__(self, data, n, w, a)
        #self.out_s1(ind=41)

        self.k = k
        self.cm = {}
        self.d = d

    def build_cm(self, t):
        
        for i in range(t):
            self.mask_space()


        self.print_cm()


        n_matrix = self.m-self.n+1
        
        print "Sparse matrix size: %d (%.2f%% reduction)"%(len(self.cm), (100.*len(self.cm))/(n_matrix*(n_matrix-1)/2))
        
        print "-----------------------------------------------"


        self.s = self.cm_threshold_s(thresh_factor=1.)*t

        print "CM Threshold s = ", self.s
        top_entr = self.sort_cm_by_value()[0]


        print "max value of cm", top_entr


        cands = self.get_motif_cand()

        
        i1, i2, v = cands[0]


        self.plot_word_comp(i1, i2, 'comp_%s_%s.png'%(i1, i2))


        #i1, i2 = 40, 26

        #self.plot_word_comp(i1, i2, 'comp_%s_%s.png'%(i1, i2))

        self.out_s1(ind=i1)
        self.out_s1(ind=i2)
        


    def mask_space(self):
        mask = random.sample(range(self.w), self.k)

        mask_dict = {}

        for col in self.s_matrix:
        #for col in range(self.m-self.n+1):
            # decompress
            # s_col = df_get_col_rle(self.s_matrix, col)

            s_col = self.s_matrix[col]

            # Triar la projeccio del window 'col'
            word = tuple(s_col[mask].tolist())
            # two or more words hash to the same bucket
            if word in mask_dict:
                self.update_cm(mask_dict[word], col)
                mask_dict[word].append(col)
            else:
                mask_dict[word] = [col]


        


    def sort_cm_by_value(self, rev=True):
        return sorted(self.cm.values(), key=lambda x: x[1],reverse=rev)
        

    def cm_threshold_s(self, thresh_factor=10.):
        # EXPECTATION VALUE
        # AL TANTO QUE ES PER ITERATION!!!!

        s = thresh_factor*expect_s_thresh(num_words=len(self.s_matrix.columns), a=len(self.symbs), w=self.w, d=self.d, mask_size=self.k)

        return s


    def update_cm(self, existing_indices, new_index):

        for old_i in existing_indices:

            key = frozenset([old_i, new_index])

            # Trivial matches work
            if np.abs(new_index - old_i) >= self.n:

                if key in self.cm:
                    self.cm[key][1] += 1
                else:
                    self.cm[key] = [key, 1]




    def get_motif_cand(self):

        sort_vals = self.sort_cm_by_value()
        

        for valors in sort_vals[:10]:
            print valors

        diagonal_dict = {}


        for v in sort_vals:
            i,j = tuple(v[0])
            diag_key = np.abs(j-i)
            # Trying to eliminate trivial matches, does it have to be strictly n?
            if diag_key > self.n:

                d1 = self.data[i:i+self.n]
                d2 = self.data[j:j+self.n]
                distancia = dist_function(d1, d2, self.d)

                if diag_key in diagonal_dict:
                    if distancia < diagonal_dict[diag_key][2]:
                        diagonal_dict[diag_key] = (i,j,distancia)

                else:
                    diagonal_dict[diag_key] = (i,j,distancia)
                
                
                

        top_pair_cand = sorted(diagonal_dict.values(), key=lambda x:x[2])


        print "Top pair", top_pair_cand[0]

        return top_pair_cand



    def plot_word_comp(self, ind1, ind2, filename='test.png'):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(211)


        d1 = self.data[ind1:ind1+self.n]
        d2 = self.data[ind2:ind2+self.n]

        ax.plot(range(ind1,ind1+self.n),d1, color='g', marker='', linestyle='-', lw=2., mew=2.)
        ax.plot(range(ind2,ind2+self.n),d2, color='b', marker='', linestyle='-', lw=2., mew=2.)
        ax.plot(range(len(self.data)), self.data, color='k', marker='', linestyle='-', alpha=.5, mew=1.)


        ax2 = fig.add_subplot(212)

        data_norm = normalize_cs(self.data)

        d1p = data_norm[ind1:ind1+self.n]
        d2p = data_norm[ind2:ind2+self.n]

        ax2.plot(range(self.n), d1p, color='g', marker='', linestyle='-', lw=2., mew=2.)
        ax2.plot(range(self.n), d2p, color='b', marker='', linestyle='-', lw=2., mew=2.)


        ax.set_xlabel('t')
        ax.set_ylabel('T')
        

        plt.suptitle('n = %s; D()=%.2f; d = %d; cm = %d'%(self.n, dist_function(d1, d2, self.d), self.d, self.cm[frozenset([ind1, ind2])][1]))


        fig.savefig(filename, format='png')
        plt.close()




    def print_cm(self, filename='cm.png'):

        from pylab import cm

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)


        x_shift = [-0.5 + e for e in range(0, self.m-self.n+1)]
        y_shift = [-0.5 + e for e in range(0, self.m-self.n+1)]

        X_shift, Y_shift = np.meshgrid(x_shift, y_shift)


        plot_array = []
        for i in range(self.m-self.n+1):
            plot_array.append([])
            
            for j in range(self.m-self.n+1):
                key = frozenset([i,j])
                if key in self.cm and i > j:
                    plot_array[i].append(self.cm[key][1])
                else:
                    plot_array[i].append(0.0)

        min_val = self.sort_cm_by_value(rev=False)[0][1]
        max_val = self.sort_cm_by_value(rev=True)[0][1]

        var = plt.pcolormesh(np.array(X_shift),np.array(Y_shift),np.array(plot_array).T, cmap=cm.gray_r, vmin=0, vmax=max_val)
        plt.colorbar(var)


        plt.xlim(-0.5, self.m-self.n+0.5)
        plt.ylim(-0.5, self.m-self.n+0.5)


        ax.set_xlabel('')
        ax.set_ylabel('')

        #plt.suptitle('Collision Matrix')
        plt.title('Collision Matrix')

        fig.savefig(filename, format='png')
        plt.close()




# k < w - d
if __name__ == '__main__':

    data = read_data('data.txt')

    params_nwa = (14,7,8)

    bar = TSP2(data, params=params_nwa, k=5, d=2)
    
    bar.build_cm(t=50)


    

































