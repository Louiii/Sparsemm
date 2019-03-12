from os import walk

class LikwidReader:
    def __init__(self, folder):
        self.folder = folder
    def make(self):
        f = []
        for (dirpath, dirnames, filenames) in walk(self.folder):
            f.extend(filenames)
            break
        f.sort()
#        f = f[8:]+f[:8]
        print(f)
        f.pop(0)
        runtimes = {}
        for fname in f:
            with open(self.folder+fname) as f_:
                content = f_.readlines()
                
                stats = {}
            for line in content:
#                print(stats)
                if 'NON-ZEROS' in line or 'NZ' in line:
                    print(line)
                    positions = [pos for pos, char in enumerate(line) if char == ':' or char == ',' or char == '=']
                    stats['Non-Zeros'] = int(line[positions[0]+1:positions[1]])
                    c = -1
                    if len(positions)>3:
                        c = positions[3]
                    print(line[positions[2]+1:c])
                    stats['n'] = int(line[positions[2]+1:c])
                if 'Runtime [s]' in line:
                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    stats['Runtime [s]'] = float(line[positions[1]+1:positions[2]])
                    if 'Runtime [s]' not in stats:
                        stats['Runtime [s]'] = float(line[positions[1]+1:positions[2]])
                    else:
                        stats['Runtime [s]'] += float(line[positions[1]+1:positions[2]])
                    
#                if 'call count' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    stats['call count'] = int(line[positions[1]+1:positions[2]])
                if 'Clock' in line:
                    positions = [pos for pos, char in enumerate(line) if char == '|']
                    stats['Clock'] = float(line[positions[1]+1:positions[2]])
                if 'CPI' in line:
                    positions = [pos for pos, char in enumerate(line) if char == '|']
                    stats['CPI'] = float(line[positions[1]+1:positions[2]])
#                if 'Load to store ratio' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    stats['Load to store ratio'] = float(line[positions[1]+1:positions[2]])
                if 'ICACHE_MISSES'  in line:
                    positions = [pos for pos, char in enumerate(line) if char == '|']
                    num = line[positions[2]+1:positions[3]]
                    if 'e' in num: 
                        n = float(num[:num.index('e')])*10**int(num[num.index('e')+1:])
                    else:
                        n = int(num)
                    if 'Cache misses' not in stats:
                        stats['Cache misses'] = int(n)
                    else:
                        stats['Cache misses'] += int(n)
#                if 'INSTR_RETIRED_ANY' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    num = line[positions[2]+1:positions[3]]
#                    if 'e' in num: 
#                        n = float(num[:num.index('e')])*10**int(num[num.index('e')+1:])
#                    else:
#                        n = int(num)
#                    stats['INSTR_RETIRED_ANY'] = ( line[positions[1]+1:positions[2]].replace(" ", ""), int(n) )
#                if 'CPU_CLK_UNHALTED_CORE' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    num = line[positions[2]+1:positions[3]]
#                    if 'e' in num: 
#                        n = float(num[:num.index('e')])*10**int(num[num.index('e')+1:])
#                    else:
#                        n = int(num)
#                    stats['CPU_CLK_UNHALTED_CORE'] = ( line[positions[1]+1:positions[2]].replace(" ", ""), int(n) )
#                if 'CPU_CLK_UNHALTED_REF' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    num = line[positions[2]+1:positions[3]]
#                    if 'e' in num: 
#                        n = float(num[:num.index('e')])*10**int(num[num.index('e')+1:])
#                    else:
#                        n = int(num)
#                    stats['CPU_CLK_UNHALTED_REF'] = ( line[positions[1]+1:positions[2]].replace(" ", ""), int(n) )
#                if 'MEM_UOPS_RETIRED_LOADS_ALL' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    num = line[positions[2]+1:positions[3]]
#                    if 'e' in num: 
#                        n = float(num[:num.index('e')])*10**int(num[num.index('e')+1:])
#                    else:
#                        n = int(num)
#                    stats['MEM_UOPS_RETIRED_LOADS_ALL'] = ( line[positions[1]+1:positions[2]].replace(" ", ""), int(n) )
#                if 'MEM_UOPS_RETIRED_STORES_ALL' in line:
#                    positions = [pos for pos, char in enumerate(line) if char == '|']
#                    num = line[positions[2]+1:positions[3]]
#                    if 'e' in num: 
#                        n = float(num[:num.index('e')])*10**int(num[num.index('e')+1:])
#                    else:
#                        n = int(num)
#                    stats['MEM_UOPS_RETIRED_STORES_ALL'] = ( line[positions[1]+1:positions[2]].replace(" ", ""), int(n) )
            
            runtimes[fname] = stats
            print(fname)
        print('COMPLETED')
        return runtimes
    
    def select(self, readings, x_var, y_var):
        l, xss, yss = [], [], []
        for matrix in readings:
            l.append(matrix)
            if type(readings[matrix][x_var]) != tuple:
                xss.append([ readings[matrix][x_var] ])
            else:
                xss.append([ readings[matrix][x_var][1] ])
            if type(readings[matrix][y_var]) != tuple:
                yss.append([ readings[matrix][y_var] ])
            else:
                yss.append([ readings[matrix][y_var][1] ])
        return xss, yss, l
        
