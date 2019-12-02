import numpy as np

"""read QE or VASP structures to generate the basis vectors and atomic positions"""

def readst(fname):
	f=open(fname)
	content = f.readlines()
	s=['','','']; scaling_factor=1.0
	if 'ATOMIC_SPECIES\n' in content:
		for i in range(len(content)):
			if 'celldm' in content[i] and not '!' in content[i]:
				scaling_factor=float(content[i].strip('\n').split('=')[1])*0.529177
			if 'CELL_PARAMETERS' in content[i]:
				s[0]=content[i+1]
				s[1]=content[i+2]
				s[2]=content[i+3]
	else:
		s[0]=content[2]
		s[1]=content[3]
		s[2]=content[4]

	for i in range(3): 
		s[i]=s[i].strip('\n')
		s[i]=s[i].split(' ')
		for j in range(s[i].count('')):
			s[i].remove('')
		s[i]=list(map((lambda x:float(x)),s[i]))
		s[i]=np.array(s[i])*scaling_factor

	alist=[]
	if 'ATOMIC_SPECIES\n' in content:
		for i in range(len(content)):
			if 'nat' in content[i] and not '!' in content[i]:
				nat=content[i]
		nat=int(str(list(filter(str.isdigit, nat))).replace(',',' ').strip('[]').replace("'",'').replace(' ',''))
	else:
		sto=content[6].strip('\n').split(' ')
		for j in range(sto.count('')):
			sto.remove('')
		sto1=[i for i in sto if i != ' ']
		sto2=[int(i) for i in sto1]
		nat=sum(sto2)
		nam=list(content[5].strip('\n').split(' '))
		nam1=[i for i in nam if i != '']
		for i in range(len(sto2)):
			for j in range(sto2[i]):
				alist.append(nam1[i])

	p=[]
	for i in range(nat):
		p.append('')

	pp=[]
	for i in range(nat):
		pp.append('')

	if 'ATOMIC_SPECIES\n' in content:
		for i in range(len(content)):
			if 'ATOMIC_POSITIONS' in content[i]:
				for j in range(nat):
					pp[j]=content[i+1+j]

		for i in range(nat): 
			pp[i]=pp[i].strip('\n')
			pp[i]=pp[i].split(' ')
			for j in range(pp[i].count('')):
				pp[i].remove('')
			alist.append(pp[i][0])
			p[i]=pp[i][1:4]
			p[i]=list(map((lambda x:float(x)),p[i]))
	else:
		for i in range(nat):
			p[i]=content[i+8]
			p[i]=p[i].strip('\n')
			p[i]=p[i].split(' ')
			for j in range(p[i].count('')):
				p[i].remove('')
			p[i]=list(map((lambda x:float(x)),p[i][0:3]))
	for k in range(len(alist)):
		alistt=''.join(i for i in alist[k] if not i.isdigit())
		alist[k]=alistt
	pp=list(zip(alist,p))
	for i in range(nat):
		pp[i]=list(pp[i])

	return nat,s,pp,p
