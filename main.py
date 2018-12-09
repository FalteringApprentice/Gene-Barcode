import numpy as np
import cv2 as cv
from itertools import product

def cal_N(k):
    if (k % 2) == 0:
        return (4**k+4**(k//2))//2
    else:
        return 4**k//2

def split_substr(s, l):
    '''
    均分序列，将序列s均分为长度为l的多个子序列
    :param s: 源序列
    :param l: 分段段长
    :return: 分段后的序列数组
    '''
    s_len=len(s)    #s的长度
    n=s_len // l    #共可得将序列分为几段
    res=[]
    ix = 0
    for i in range(n):
        res.append(s[ix:ix+l])
        ix += l
    return res


def split_sublist(list,L):
    '''
    分割数组list为L个子数组
    :param list:待分割数组
    :param L:分割份数
    :return:分割后的数组
    '''
    res=[]
    list_len=len(list)
    piece_len=list_len//L  #分割每份长度
    left_len=list_len-L*piece_len
    ix=0
    for i in range(L):
        if i<left_len:
            res.append(list[ix:ix+piece_len+1])
            ix+=piece_len+1
        else:
            res.append(list[ix:ix + piece_len])
            ix += piece_len
    return res

def reverse_complement(seq):
    '''
    求反向互补序列
    :param seq: 原序列
    :return: 反向互补序列
    '''
    comp={'A':'T','T':'A','C':'G','G':'C'}
    res=''
    for s in seq:
        res+=comp[s]
    return res[::-1]


def generate_kmer(k):
    '''
    生成kmer核苷酸排列数组
    :param k: 参数k
    :return: kmer排列数组
    '''
    pairs=[]
    kmers = [''.join(i) for i in product('AGTC', repeat=k)]
    for i in kmers:
        rc=reverse_complement(i)
        for j in kmers:
            if j==rc and j!=rc[::-1]:
                kmers.remove(j)
                pairs.append((i,rc))
            elif j==rc and j==rc[::-1]:
                pairs.append(i)
    return pairs


def read():
    '''
    数据读取
    :param seq:存入的字典
    :return:无
    '''
    seqs={}
    with open('hs_ref_GRCh38.p7_chr1.fa') as f:
        curName=''
        curSeq=''
        buffer=f.read()
        buffer=buffer.split('>')
        buffer.pop(0)
        for seq in buffer:
            bodybegin=seq.find('\n')
            meta=seq[0:bodybegin]
            curName=meta.split()[0].strip()
            curSeq=seq[bodybegin+1:-1].replace('\n','').strip()
            seqs[curName]=curSeq
        f.close()
    return seqs


def cal_freq(seq,s):
    '''
    计算kmer在seq中出现的频率
    :param seq:长序列
    :param s:短序列字符串或短序列字符串与它的反向互补序列构成的元组
    :return:s在seq中出现的频率
    '''
    cnt=0
    l=len(seq)
    s_len=0
    if isinstance(s,str):
        s_len=len(s)
        for i in range(len(seq)-s_len):
            if seq[i:i+s_len]==s:
                cnt=cnt+1
    elif isinstance(s,tuple):
        s_len=len(s[0])
        for i in range(len(seq)-s_len):
            if seq[i:i+s_len]==s[0] or seq[i:i+s_len]==s[1]:
                cnt=cnt+1
    return cnt/(len(seq)-s_len+1)


def map_freq_to_gray(freq,k,L=14):
    '''
    将频率映射为灰度阶
    :param freq: 频率矩阵
    :param k: k-mer
    :param L: 分段长度
    :return: 值对应的灰度字典
    '''
    N=cal_N(k)
    rownum=len(freq)
    res=[]
    for i in range(rownum):
        cur_map={}
        cur_row=freq[i]
        S=np.sort(cur_row)
        avg_S=np.average(S)
        splited_list=split_sublist(S,L)
        gradient=256/L
        for j in range(splited_list.__len__()):
            sub_list=splited_list[j]
            list_map={}
            for k in range(sub_list.__len__()):
                list_map[sub_list[k]]=int(j*gradient)
            for item in list_map:
                if (not cur_map.get(item)) and item!=0.0:
                    cur_map[item]=list_map[item]
                elif (item==0.0):
                    cur_map[item]=0
        res.append(cur_map)
    return res


def barcode(seq,M=1000,k=4):
    '''
    生成基因条形码
    :param seq: 序列
    :param M: 分段段长
    :param k: k-mer的参数k
    :return:无
    '''
    N=cal_N(k)
    rownum = seq.__len__() // M
    colnum = N
    freq=[]
    bps=split_substr(seq,M)
    kmers = generate_kmer(k)
    for i in range(rownum):
        row=[]
        for j in range(colnum):
            cur_seq=bps[i]
            row.append(cal_freq(cur_seq,kmers[j]))
        freq.append(row)
    gray_map=map_freq_to_gray(freq.copy(),k)

    img = np.zeros([rownum, colnum, 1], dtype=np.uint8)
    for i in range(rownum):
        for j in range(colnum):
            img[i][j][0] = gray_map[i][freq[i][j]]
    cv.imshow("gene barcode",img)
    cv.waitKey(0)
    cv.destroyAllWindows()


def main():
    seqs = read()
    for id in seqs:
        barcode(seqs[id])


if __name__ =='__main__':
    main()
