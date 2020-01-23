import sys
import random
import math
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
import numpy as np

# execution style : python < name of file > < number of frames> < number of elements in each frame >

def ReadData(filename,a,NOF):
    # Read the file, splitting by lines
    
    

    y = 0
    f = open(filename, "r")
    lines = f.readlines()

    result = []
    for x in lines:
        result.append(x)
        y = y + 1

    f.close()


    counter = 0
    refine = []
    x = 0

    while x <= a * (NOF + 1):
        if x == 0 or x == 1:  # skip first 2 lines
            counter = counter + 1
            x = x + 1
            # print("x is %s and counter in %s in stage 1" % (x, counter))

        elif 2 <= counter <= a + 1:
            refine.append(result[x])
            x = x + 1
            counter = counter + 1
            # print("x is %s and counter in %s in stage 2" % (x, counter))

        elif counter == a + 2:
            counter = 2
            x = x + 3
            if x > y:
                break
        # print("x is %s and counter in %s in stage 3" % (x, counter))

        else:
            print(0)

    file_out = open("array.txt", "w")

    # admit values into a large 2D array

    rows, cols = (a * NOF, 4)
    array = [[0 for i in range(cols)] for j in range(rows)]

    for i in range(a * NOF):
        t = refine[i][0:8]
        id = refine[i][15:20]  # Column values may vary on the basis of gromacs version
        x = refine[i][21:28]
        y = refine[i][29:36]
        z = refine[i][37:44]
        array[i][0] = id
        array[i][1] = x
        array[i][2] = y
        array[i][3] = z

    for i in range(rows):
        file_out.writelines(str(array[i]) + '\n')
        # view array.txt to check
    file_out.close()

    # 2d big to 3D small conversion

    rows, cols, pages = (a, 4, NOF)
    array3d = [[[0 for k in range(pages)] for i in range(cols)] for j in range(rows)]

    h = 0

    for k in range(pages):
        for i in range(rows):
            array3d[i][0][k] = float(array[h][0])
            array3d[i][1][k] = float(array[h][1])
            array3d[i][2][k] = float(array[h][2])
            array3d[i][3][k] = float(array[h][3])
            h = h + 1

    file_out1 = open("array3d.txt", "w")

    # view array3d.txt to check

    for k in range(pages):
        for i in range(rows):
            file_out1.writelines(
                str(array3d[i][1][k]) + ',' + str(array3d[i][2][k]) + ',' + str(array3d[i][3][k]) + '\n')

    file_out1.close()

    fileName = "array3d.txt"
    f = open(fileName, 'r')
    lines = f.read().splitlines()
    f.close()

    items = []

    for i in range(1, len(lines)):
        line = lines[i].split(',')
        itemFeatures = []

        for j in range(len(line)):
            v = float(line[j])  # Convert feature value to float
            itemFeatures.append(v)  # Add feature value to dict

        items.append(itemFeatures)

    random.shuffle(items)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    return items;


def FindColMinMax(items):
    n = len(items[0]);

    minima = [sys.maxsize for i in range(n)];
    maxima = [-sys.maxsize - 1 for i in range(n)];

    for item in items:
        for f in range(len(item)):

            if (item[f] < minima[f]):
                minima[f] = item[f];

            if (item[f] > maxima[f]):
                maxima[f] = item[f];

    # print(maxima, minima)
    return minima, maxima;


def InitializeMeans(items, k, cMin, cMax):
    # Initialize means to random numbers between
    # the min and max of each column/feature
    f = len(items[0]);  # number of features
    means = [[0 for i in range(f)] for j in range(k)];

    for mean in means:
        for i in range(len(mean)):
            # Set value to a random float
            # (adding +-1 to avoid a wide placement of a mean)
            mean[i] = random.uniform(cMin[i] + 1, cMax[i] - 1);

    # print(means)
    return means;


def EuclideanDistance(x, y):
    S = 0  # The sum of the squared differences of the elements
    for i in range(len(x)):
        S += math.pow(x[i] - y[i], 2)

    return math.sqrt(S)  # The square root of the sum


def UpdateMean(n, mean, item):
    for i in range(len(mean)):
        m = mean[i];
        m = (m * (n - 1) + item[i]) / float(n);
        mean[i] = round(m, 3);

    return mean;


def Classify(means, item):
    # Classify item to the mean with minimum distance
    minimum = sys.maxsize;
    index = -1;

    for i in range(len(means)):

        # Find distance from item to mean
        dis = EuclideanDistance(item, means[i]);

        if (dis < minimum):
            minimum = dis;
            index = i;

    return index;


def CalculateMeans(k, items, maxIterations=100000):
    # Find the minima and maxima for columns
    cMin, cMax = FindColMinMax(items);

    # Initialize means at random points
    means = InitializeMeans(items, k, cMin, cMax);

    # print('means \n')
    # print(means)

    # Initialize clusters, the array to hold
    # the number of items in a class
    clusterSizes = [0 for i in range(len(means))];

    # An array to hold the cluster an item is in
    belongsTo = [0 for i in range(len(items))];

    # Calculate means
    for e in range(maxIterations):

        # If no change of cluster occurs, halt
        noChange = True;
        for i in range(len(items)):

            item = items[i];

            # Classify item into a cluster and update the
            # corresponding means.
            index = Classify(means, item);

            clusterSizes[index] += 1;
            cSize = clusterSizes[index];
            means[index] = UpdateMean(cSize, means[index], item);

            # Item changed cluster
            if (index != belongsTo[i]):
                noChange = False;

            belongsTo[i] = index;

        # Nothing changed, return
        if (noChange):
            break;

    return means;


def FindClusters(means, items):
    clusters = [[] for i in range(len(means))];  # Init clusters

    for item in items:
        # Classify item into a cluster
        index = Classify(means, item);

        # Add item to cluster
        clusters[index].append(item);

    return clusters;


def main():


    ###### command line section ####################

    filename = sys.argv[1]
    print('filename is %s' % (filename))
    NOF = int(sys.argv[2])
    print('NOF is %s' % (NOF))
    a = int(sys.argv[3])
    print('No. of elements in a single frame is %s' % (a))

    #####################################################################
    store_dia = [0 for r in range(3)]
    for count in range(3):


        np.random.seed(19680801)
        print('Doing K Means of 2 in experiment %s' % (count+1))
        items = ReadData(filename,a,NOF)
        # print(items)

        means = CalculateMeans(2, items)

        clusters = FindClusters(means, items)

        centroid_x = [0 for r in range(2)]
        centroid_y = [0 for r in range(2)]
        centroid_z = [0 for r in range(2)]

        for i in range(len(means)):
            centroid_x[i] = float(means[i][0])
            centroid_y[i] = float(means[i][1])
            centroid_z[i] = float(means[i][2])

        ##################################################
        np.random.seed(19680801)
        print('Doing K Means of 1 to find COM in experiment %s '% (count+1))
        items = ReadData(filename,a,NOF)
        # print(items)

        means = CalculateMeans(1, items)

        clusters = FindClusters(means, items)

        com_x = [0 for r in range(2)]
        com_y = [0 for r in range(2)]
        com_z = [0 for r in range(2)]

        for i in range(len(means)):
            com_x = float(means[i][0])
            com_y = float(means[i][1])
            com_z = float(means[i][2])

        diameter_1 = math.sqrt(
            math.pow(centroid_x[0] - com_x, 2) + math.pow(centroid_y[0] - com_y, 2) + math.pow(centroid_z[0] - com_z,
                                                                                               2))
        diameter_2 = math.sqrt(
            math.pow(centroid_x[1] - com_x, 2) + math.pow(centroid_y[1] - com_y, 2) + math.pow(centroid_z[1] - com_z,
                                                                                               2))
        store_dia[count] = (diameter_1 + diameter_2)
        print("Result : % s \n" % (store_dia[count]))



    sum_of_dia = ( store_dia[0] + store_dia[1] + store_dia[2] )
    avg_diameter=sum_of_dia/3
    print(avg_diameter)

    file_out1 = open("check_dump.txt", "w")

    # check_dump.txt to check


    file_out1.write(str(avg_diameter))

    file_out1.close()
    
    print("check output in check_dump.txt \n")


if __name__ == "__main__":
    main()
