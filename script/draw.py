s = 's,s,s,s'
s.strip

def lineToXYZ(line):
    '''
    Input: a line from df.readline()
    Output: a list of [x,y,z] float value
    '''
    line = line.strip("\n")
    line = line.split(",")
    for i in range(3):
        line[i] = float(line[i])
    return line


file_name = './pos.txt'
df = file(file_name, 'r')
print lineToXYZ(df.readline())

def main():
    '''
    need 1+3+3+3 vectors
    t, pos, vel and acc
    '''
    time = []
    pos_x = []
    pos_y = []
    pos_z = []
    vel_x = []
    vel_y = []
    vel_z = []
    acc_x = []
    acc_y = []
    acc_z = []

    '''
    open 4 files, time, pos, vel, acc
    '''
    timef = file("./time.txt",'r')
    posf = file("./pos.txt",'r')
    velf = file("./vel.txt",'r')
    accf = file("./acc.txt",'r')

    'read the file until end'
    time = timef.readline()
    while not time == '':
        time = float(time.strip('\n')[0])
        pos = lineToXYZ(posf.readline())
        vel = lineToXYZ(velf.readline())
        acc = lineToXYZ(accf.readline())
        print time, pos, vel, acc

    


if __name__ == '__main__':
    main()