'''
A set of functions that computes the entropy, either shannon or kolmogorov, of
a string. The kolmogorov complexity is an approximation using zlip and shannon
was taken off of activestate and modified to use variable wordsizes. 

@author: dbgoodman
'''
## {{{ http://code.activestate.com/recipes/577476/ (r1)
# Shannon Entropy of a string
# = minimum average number of bits per symbol
# required for encoding the string
#
# So the theoretical limit for data compression:
# Shannon Entropy of the string * string length
# FB - 201011291
import math
import sys
import zlib
import warnings
import random
import itertools

def shannon(st, wordsize=1):
    #st = 'aabcddddefffg' # input string
    # Shannon entropy for this would be 1 bit/symbol
    #st = '00010101011110' 
    
    #print 'Input string:'
    #print st
    #print 
    
    wordsize = int(wordsize)
    #print 'Word size:'
    #print wordsize
    #print 
    
    
    if wordsize == 1:
        stList = list(st)
    else:
        stList = \
            [st[i:i+wordsize] for i in range(0, len(st), wordsize)]
    
    alphabet = list(set(stList)) # list of symbols in the string
    
    #print 'Alphabet of symbols in the string:'
    #print alphabet
    #print
    
    # calculate the frequency of each symbol in the string
    freqList = []
    for symbol in alphabet:
        ctr = 0
        for sym in stList:
            if sym == symbol:
                ctr += 1
        freqList.append(float(ctr) / len(stList))
    #print 'Frequencies of alphabet symbols:'
    #print freqList
    #print
    # Shannon entropy
    ent = 0.0
    for freq in freqList:
        ent = ent + freq * math.log(freq, 2)
    ent = -ent
    return ent
    #print 'Shannon entropy:'
    #print ent
    #print 'Minimum number of bits required to encode each symbol:'
    #print int(math.ceil(ent))
## end of http://code.activestate.com/recipes/577476/ }}}

def kolmogorov(st):
    '''
    use zlib.compress to approximate the kolmogorov score. the difference be-
    -tween the compressed and original is 8 bytes, so we subtract that from the
    compressed length. this approximation isn't valid for strings less than 5.
    As it approaches 1, it means the string is incompressible (high entropy. As
    it approaches 0, it means the string has very low entropy.
    '''
    if len(st) < 5:
        warnings.warn('''Kolmogorov approximation is not valid for strings 
                         smaller than len 5''')
    l = float(len(st))
    compr = zlib.compress(st)
    c = float(len(compr)-8)
    return c/l

def test_kolmogorov(seq_size= 10, seqs= 100, k='['+'0.25, '*4+']'):
    
    k = eval(k)
    k = map(int, k)
    
    seq_size = int(seq_size)
    seqs = int(seqs)
    
    sample_space = ''.join(map(lambda a: a[0]*a[1], zip('ACGT',k)))
    
    rand_sum = 0
    repeat_sum = 0
    
    for i in range(seqs):
        rand_seq = ''.join([random.choice('AGTC') 
                                    for x in range(seq_size)])
        
        rand_kol = kolmogorov(rand_seq)
        
        #print 'Rand:\t{0}\t{1}'.format(rand_seq, rand_kol)
        
        rand_sum = rand_sum + rand_kol
        
        repeat_seq = ''.join([random.choice(sample_space) 
                                    for x in range(seq_size)])
        
        repeat_kol = kolmogorov(repeat_seq[0:seq_size])
        
        #print 'Rep:\t{0}\t{1}'.format(repeat_seq[0:seq_size], repeat_kol)
        
        repeat_sum = repeat_sum + repeat_kol
                       
    rand_sum = rand_sum / seqs
    repeat_sum = repeat_sum / seqs
    
    #print 'avg kolmogorov entropy of random sequences: %f' % rand_sum
    #print 'avg kolmogorov entropy of repeat sequences: %f' % repeat_sum
    return (rand_sum,repeat_sum)

def test_kolmogorov_xy(start,end,step,k):
    x = range(*map(int,[start,end,step]))
    for length in x:
        rand_sum,repeat_sum = test_kolmogorov(length,1000,k)
        print '\t'.join(map(str,[length, rand_sum, repeat_sum]))


if __name__ == "__main__":
    if sys.argv[1] in locals():
        out = eval(sys.argv[1]+'(*sys.argv[2:])')
        print out
    else:
        raise(ValueError('no function called %s' %sys.argv[1]))
