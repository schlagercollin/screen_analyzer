import csv
import openpyxl

wb = openpyxl.Workbook()
ws = wb.active
ws.title = "Library Data"

f = open('/Users/collinschlager/Documents/Rohatgi_Lab/screen_analyzer/tmp/data/output/debug_large_matched.csv')
reader = csv.reader(f, delimiter=',')
for row in reader:
    ws.append(row)
f.close()

wb.save('test.xlsx')
