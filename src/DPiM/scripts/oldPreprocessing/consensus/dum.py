agg = aggregator() # numpy happens inside this class
for f in oldJsonFiles:
    sparseMat = readJson(f) # json parsing happens here
    agg.push(obj.sparseMat) 

oldMeans = agg.means() 

for f in newJsonFiles:
    sparseMat = readJson(f)
    agg.push(sparseMat)

newMeans = agg.means()
stddev = agg.stddev()
maxPair = findMax(stdDev) # find the indices (protein pair) with the highest sd

if (abs(newMeans.at(maxPair) - oldMeans.at(maxPair))/stddev.at(maxPair)
    < limit):
    return success
else:
    return fail

# This class calculates mean and variance on a running basis
# It's essentially just simple arithmetic
# I'm implementing this as though what's added isn't a matrix but a single value
class aggregator(object):
    def __init__(self):
        self.count=0
        ## the following should be (floating point) matrices
        self.prevMean=0 
        self.newMean=0
        self.prevSqr=0
        self.newSqr=0

    # add a new value, update the statistics
    def push(self, val):
        self.count = self.count+1
        c = self.count
        if (c == 1):
            self.prevMean=self.newMean=val
            self.oldSqr = 0
        else:
            self.newMean = self.prevMean + (val-self.oldMean)/c
            self.newSqr = self.prevSqr + \
                          (val - self.oldMean)*(val - self.newMean) 
            self.oldMean=self.newMean
            self.oldSqr=self.oldSqr

    def means(self):
        # check if self.count == 0
        return self.newMean

    def stddev(self):
        # check if self.count == 0
        var = self.newSqr/(self.count-1)
        return sqrt(var)
