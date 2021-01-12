


## some params related to exits
NEXITS = 5 # how many cashiers in the store
CASHIERD = 14 # distance between cashiers
ENTRANCEPOS = 0 # x-coord of entrance to the store
EXITPOS = 1 # Lx - y-coord of the first cashier 

## diffusion coeffs
DIFFCOEFF = 5e-2
ACSINKCOEFF= 1e-2 # 1 / tau
PLUMEMIN=0.0

## plume spreading parameters
PLUMELIFETIME = 20 ## lifetime of plume for discrete plumes without diffusion
PLUMECONCINC = 40000.0 ## aerosol concentration in coughing event
PLUMECONCCONT = 5.0 ## continuous aerosol emission
PLUMEPLOTMIN=1 ## paramter for plotting method

CASHIERTIMEPERITEM = 1 ## waiting time multiplier on cashier
BLOCKRANDOMSTEP=0.8  ## parameter giving prob that customer takes random step if path blocked (i.e. another customer blocking)
PROBSPREADPLUME =  1./3600 # i.e. prob of cough per sec

EXPOSURELIMIT=1 ## a threshold for counting advanced exposure statistics 
## limits for maximum waiting time when a target from the shopping list is found
MAXWAITINGTIME=2 
MINWAITINGTIME=1

## some simulation parameters for the customer behaviour
MAXSHOPPINGLIST=20 
WEIRDQUEUELIMIT = 39 ## parameter for starting plotting images when queues grow larger than the value


