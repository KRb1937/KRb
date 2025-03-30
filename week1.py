ref_str="CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"
que_str="CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"


origin_base="ATCG"
reverse_base="TAGC"
trans_base=str.maketrans(origin_base,reverse_base)

rev_str=que_str.translate(trans_base)
rev_str=rev_str[::-1]

def do_hash(str,mod):
    prehash=[0 for i in range(len(str)+1)]
    for i in range(len(str)):
        prehash[i+1]=(prehash[i]*257+ord(str[i]))%(mod)
    return prehash

def do_power(str,mod):
    prepower=[1 for i in range(len(str)+1)]
    for i in range(len(str)):
        prepower[i+1]=(prepower[i]*257)%(mod)
    return prepower

ref_prehash1=do_hash(ref_str,1000000007)
que_prehash1=do_hash(que_str,1000000007)
rev_prehash1=do_hash(rev_str,1000000007)
ref_prehash2=do_hash(ref_str,1000000009)
que_prehash2=do_hash(que_str,1000000009)
rev_prehash2=do_hash(rev_str,1000000009)

ref_prepower1=do_power(ref_str,1000000007)
que_prepower1=do_power(que_str,1000000007)
rev_prepower1=do_power(rev_str,1000000007)
ref_prepower2=do_power(ref_str,1000000009)
que_prepower2=do_power(que_str,1000000009)
rev_prepower2=do_power(rev_str,1000000009)

def get_hash(prehash,prepower,begin,end,mod):
    return (prehash[end]-prehash[begin]*prepower[end-begin]%mod)%mod



ref_hash_dict={}
for i in range(1,len(ref_str)+1):
    for j in range(len(ref_str)-i+1):
        hash_value=(get_hash(ref_prehash1,ref_prepower1,j,i+j,1000000007),get_hash(ref_prehash2,ref_prepower2,j,i+j,1000000009))
        if hash_value not in ref_hash_dict:
            ref_hash_dict[hash_value] = []
        ref_hash_dict[hash_value].append(j)


i=0
repeat=[]

while i<len(que_str):
    repeat_len=0
    repeat_pos=0
    reverse_len=0
    reverse_pos=0
    jtop=min(len(ref_str),len(que_str)-i)
    jbot=0
    while jtop-jbot>1:
        j=(jtop+jbot)//2
        que_hash1=get_hash(que_prehash1,que_prepower1,i,i+j,1000000007)
        que_hash2=get_hash(que_prehash2,que_prepower2,i,i+j,1000000009)
        if (que_hash1,que_hash2) in ref_hash_dict:
            repeat_len=j
            repeat_pos=ref_hash_dict[(que_hash1,que_hash2)][0]
            jbot=j
        else:
            jtop=j
    jtop=min(len(ref_str),len(que_str)-i)
    jbot=0
    while jtop-jbot>1:
        j=(jtop+jbot)//2
        rev_hash1=get_hash(rev_prehash1,rev_prepower1,len(que_str)-i-j,len(que_str)-i,1000000007)
        rev_hash2=get_hash(rev_prehash2,rev_prepower2,len(que_str)-i-j,len(que_str)-i,1000000009)
        if (rev_hash1,rev_hash2) in ref_hash_dict:
            reverse_len=j
            reverse_pos=ref_hash_dict[(rev_hash1,rev_hash2)][0]
            jbot=j
        else:
            jtop=j
    larger=max(repeat_len,reverse_len)
    if larger>0:
        if repeat_len>reverse_len:
            end_pos=repeat_pos+repeat_len
            repeat.append((end_pos,repeat_len,"False"))
            i+=repeat_len
        else:
            end_pos=reverse_len+reverse_pos
            repeat.append((end_pos,reverse_len,"True"))
            i+=reverse_len
    else:
        i+=1

result={}
ref_max=0
for data in repeat:
    if (data[0],data[1]) not in result:
        if data[0]<=ref_max:
            result[(data[0],data[1])]=[data[0],data[1],1,data[2]]
        else:
            ref_max=data[0]
    else:
        result[(data[0],data[1])][2]+=1

print("|  index  |  POS in REF  |  Repeat size  |  Repeat count  |  Inverse  |")
i=1
for data in result.values():
    print("{:^1}{:^9}{:^1}{:^14}{:^1}{:^15}{:^1}{:^16}{:^1}{:^11}{:^1}".format("|",i,"|",data[0],"|",data[1],"|",data[2],"|",data[3],"|"))
    i+=1