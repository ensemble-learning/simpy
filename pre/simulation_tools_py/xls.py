PRO = {}
f = open('pro.csv', 'r')
for i in f:
    tokens = i.split(',')
    PRO[tokens[0]] = []
    PRO[tokens[0]].append(tokens[1])
    PRO[tokens[0]].append(tokens[2].strip())
f.close()

from xlwt import *
book = Workbook()
sheet1 = book.add_sheet(u'Al_nmr')
head = ['fragment', 'mol_structure', 'charges', 'pos', 'shielding parameter']

al = Alignment()
al.horz = Alignment.HORZ_CENTER
al.vert = Alignment.VERT_CENTER

borders = Borders()
borders.bottom = Borders.THICK

style = XFStyle()
style.alignment = al
style.borders = borders
row0 = sheet1.row(0)

for i, text in enumerate(head):
    row0.write(i, text, style = style)

f = open('pnmr.csv', 'r')
filename = ''
counter = 1
for i in f:
    rown = sheet1.row(counter)
    tokens = i.split(',')
    if filename == tokens[0]:
        rown.write(0, '')
    else:
        filename = tokens[0]
        molname = filename.split(r'/')[-1][:-5]
        rown.write(0, molname)
        rown.write(1, PRO[molname][1])
        rown.write(2, int(PRO[molname][0]))
    rown.write(3, int(tokens[1]))
    rown.write(4, float(tokens[2]))
        
    counter += 1
        
f.close()

book.save('tmp.xls')