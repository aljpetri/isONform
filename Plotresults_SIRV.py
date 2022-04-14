import matplotlib.pyplot as plt
import pylab
import numpy as np
from statistics import mean
def plot_scatter(nr_nodes,list):
    import matplotlib.lines as mlines
    ax = plt.gca()
    #for xe, ye in zip(getkeysList(nr_nodes), getValsList(nr_nodes)):
    #    ax.scatter([xe] * len(ye), ye)
    #x = np.linspace(*ax.get_xlim())
    #ax.plot(x, x)
    ax.yaxis.set_major_locator(plt.MaxNLocator(15))
    #plt.plot(list, list, c='red')
    plt.figtext(.8, .9, "Diversity: 20%")
    plt.title("Isoforms detected over actual Isoforms \n preprocessing isONcorrect")
    ax.scatter(getkeysList(nr_nodes),getValsList(nr_nodes), c='blue', edgecolors='none')
    #ax.yaxis.set_major_locator(plt.MaxNLocator(15))
    #set log scales for both axes
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #label the plot to make the figure better comprehensible
    ax.set_xlabel("# Isoforms Simulated")
    ax.set_ylabel("# Isoforms found by IsONform (5 experiments)")

    pylab.show()
def plot_box(resultsdict):
    labels, data = resultsdict.keys(), resultsdict.values()
    plt.figtext(.8, .9, "Diversity: 20%")
    plt.title("Isoforms detected over actual Isoforms \n preprocessing isONcorrect")
    plt.boxplot(data)
    plt.xticks(range(1, len(labels) + 1), labels)
    ax = plt.gca()
    ax.set_xlabel("# Isoforms Simulated")
    ax.set_ylabel("# Isoforms found by IsONform ( 5 experiments)")
    plt.show()
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
    filename="resultserror1.tsv"
    data=load_tsv(filename)
    resultsdict={}
    updatedresultsdict={}
    for entry in data:
        line=entry.split("\t")
        #line=line.replace("\n","")
        print(line)
        indicator=line[0].replace(" ", "")
        if indicator.isdigit():
            firstentry=int(line[0])
            print(firstentry)
            if not firstentry in resultsdict:
                values=[]
                print(line[1])
                if not int(line[1])==-1:
                    values.append(int(line[1]))
                    resultsdict[firstentry]=values
            else:
                values=resultsdict[firstentry]
                print(line[1])
                if not int(line[1]) ==-1:
                    values.append(int(line[1]))
                    resultsdict[firstentry] = values
        else:
            continue
    print("resdict",resultsdict)
    keys=resultsdict.keys()
    max_key = max(keys)
    print("max_key",max_key)
    l = list(range(2, max_key+1))
    for actualnumber, listofresults in resultsdict.items():
        #listofresults = [int(i) for i in listofresults]
        #filtered = list(filter(lambda result: result != -1, listofresults))
        print("List")
        print(listofresults)
        #print(filtered)
        avg=mean(listofresults)
        print("avg for "+str(actualnumber)+": "+str(avg))
        updatedresultsdict[actualnumber]=avg
    plot_scatter(updatedresultsdict,l)
    plot_box(resultsdict)
    # or backwards compatable


if __name__ == '__main__':
        main()
