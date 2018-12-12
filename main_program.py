import os
import struct
import imageio
import matplotlib.pyplot as plt

def chConverter(file_):
    file_.seek(0x127C)
    del_ab = struct.unpack('>d', file_.read(8))[0]

    data = []

    file_.seek(0x1800)
    while True:
        x, nrecs = struct.unpack('>BB', file_.read(2))
        # print("%a,\t%s,\t%s"%(file_.tell(),x,nrecs))
        if x == 0 and nrecs == 0:
            break
        for _ in range(nrecs):
            inp = struct.unpack('>h', file_.read(2))[0]
            if inp == -32768:
                inp = struct.unpack('>i', file_.read(4))[0]
                data.append(del_ab * inp)
            elif len(data) == 0:
                data.append(del_ab * inp)
            else:
                data.append(data[-1] + del_ab * inp)

    return data

    # with open("output_file.csv","w") as tempfile:
    #     for point in data:
    #         tempfile.write(str(point) + "\n")

def findPeaks(xdata, ydata, peakwidth=300):
    # TODO: iplement collections.deque() class for window shifting...
    # returns a list of tuples with each tuple containing (x, y) values
    peakpoints = []
    max_iters = 0
    inc_run = 0
    local_max = (xdata[0], ydata[0])
    for i in range(1, len(ydata)):
        if ydata[i] > local_max[1]:
            local_max = (xdata[i], ydata[i])
            max_iters = 0
            inc_run = 0
        elif ydata[i] <= local_max[1]:
            max_iters += 1
        if ydata[i] > ydata[i - 1]:
            inc_run += 1

        if max_iters > peakwidth and inc_run > peakwidth:
            if not (local_max[1] < 0):
                peakpoints[0].append(local_max[0])
                peakpoints[1].append(local_max[1])
            local_max = (xdata[i], ydata[i])
            max_iters = 0
            inc_run = 0
    return peakpoints


def makeGraph(pic_name, xdata, ydata, title, xlab, ylab):
    # print(xdata)

    plt.plot(xdata, ydata, "black")
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    #plt.ylim(-2, max(ydata))
    plt.ylim(-2, 100)

    plt.xlim(8, 25)
    plt.savefig(pic_name)
    # plt.show()
    plt.clf()

def analyzeData(input_directory, output_directory, photo_directory, channel_to_analyze):

    for directory in [output_directory, photo_directory]:
        for file in os.listdir(directory):
            os.remove(directory + "/" + file)

    data_files = os.listdir(input_directory)
    names = [""]
    lines = []
    float_data = []
    image_list = []

    for j, run in enumerate(data_files):
        names.append(run[11:-2])
        with open("%s/%s.csv"%(output_directory, run[11:-2]), "w") as ofile:
            with open("%s/%s/%s"%(input_directory, run, channel_to_analyze), "rb") as ifile:
                point_list = chConverter(ifile)
                for i, point in enumerate(point_list):
                    ofile.write("%s,%s\n"%(str(i/150.), str(point)))
                    try:
                        lines[i].append(str(point))
                        float_data[i].append(point)
                    except:
                        lines.append([str(len(lines)/150.), str(point)])
                        float_data.append([len(lines)/150., point])

                xdata = [time/150 for time in range(len(point_list))]
                ydata = point_list
                image_path = "%s/%s.png"%(photo_directory, run[11:-2])
                makeGraph(image_path, xdata, ydata, run[11:-2], "Time (min)", "A[280nm] (mAu)")
                image_list.append(imageio.imread(image_path))
                image_list.append(imageio.imread(image_path))
        print("%s run is converted"%run)

    imageio.mimsave("SEC-Gif.gif", image_list)
    print("Gif is made\n")

    with open("Samples Data.csv", "w") as ofile:
        ofile.write(",".join(names) + "\n")
        for line in lines:
            ofile.write(",".join(line) + "\n")

def main():
    input_directory = "data_files"
    output_directory = "output_files"
    photo_directory = "photos"
    channel_to_analyze = "DAD1A.ch"
    analyzeData(input_directory, output_directory, photo_directory, channel_to_analyze)

if __name__ == '__main__':
    main()