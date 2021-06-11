import matplotlib.pyplot as plt
import pylab
import numpy as np
from statistics import mean
def plot_scatter(nr_nodes):

    ax = plt.gca()
    #for xe, ye in zip(getkeysList(nr_nodes), getValsList(nr_nodes)):
    #    ax.scatter([xe] * len(ye), ye)
    #x = np.linspace(*ax.get_xlim())
    #ax.plot(x, x)
    #ax.yaxis.set_major_locator(plt.MaxNLocator(15))
    ax.scatter(getkeysList(nr_nodes),getValsList(nr_nodes), c='blue', edgecolors='none')
    ax.yaxis.set_major_locator(plt.MaxNLocator(15))
    #set log scales for both axes
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #label the plot to make the figure better comprehensible
    ax.set_xlabel("# Isoforms")
    ax.set_ylabel("# Isoforms found by IsONform(avg 5 tries)")

    pylab.show()
def load_tsv(filename):
    file1 = open(filename, 'r')
    elems = file1.readlines()#[1:]
    newelems=list(filter(('-e\\n\n').__ne__, elems))
    print(newelems)
    return newelems


"Helper functions to get lists of keys and values from a dict"
def getkeysList(dict):
    return dict.keys()
def getValsList(dict):
    return dict.values()

def main():
    filename="resultserror7.tsv"
    data=load_tsv(filename)
    resultsdict={}
    updatedresultsdict={}
    for entry in data:
        line=entry.split("\t")
        #line=line.replace("\n","")
        print(line)
        firstentry=line[0]
        if not firstentry in resultsdict:
            values=[]
            values.append(line[1])
            resultsdict[firstentry]=values
        else:
            values=resultsdict[firstentry]
            values.append(line[1])
            resultsdict[firstentry] = values
    for actualnumber, listofresults in resultsdict.items():
        listofresults = [int(i) for i in listofresults]
        print("List")
        print(listofresults)
        avg=mean(listofresults)
        print("avg for "+str(actualnumber)+": "+str(avg))
        updatedresultsdict[actualnumber]=avg
    plot_scatter(updatedresultsdict)


if __name__ == '__main__':
        main()
