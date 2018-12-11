def read():
    '''
    数据读取
    :return:读取的序列字典
    '''
    seqs={}
    with open('hs_ref_GRCh38.p7_chr1.fa') as f:
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

def main():
    seqs=read()
    for seq in seqs:
        out=open(seq.split('|')[-2],'w')
        out.write(seqs[seq])
        out.close()

if __name__ =='__main__':
    main()

