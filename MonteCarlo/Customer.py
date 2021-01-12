import numpy as np
from Params import *

DIRECTIONS = [np.array([-1,0]), np.array([0,1]),np.array([1,0]),np.array([0,-1])]

class Customer:


	def __init__(self, x,y, infected=0, probSpreadPlume=PROBSPREADPLUME):
		self.x =x ## initial position
		self.y =y  
		self.infected=infected ## int boolean
		self.shoppingList=[] ## list of points to visit
		self.path=None
		self.probSpreadPlume = probSpreadPlume
		self.exposure=0
		self.exposureTime=0
		self.exposureTimeThres=0
		self.timeInStore = 0
		self.initItemsList = None
		self.cashierWaitingTime = None
		self.waitingTime = 0
		self.headingForExit = 0

	## adds coordinate to the shopping list
	def addTarget(self, target):
		self.shoppingList.append(target)
		return
		
	#~ The costumer is assumed to behave logically, i.e. he/she always advances towards the closest target, which is updated utilizing this function
	def updateFirstTarget(self,store):
		shortestDist = 1.0e8
		shortInd = None
		startInd = store.getIndexFromCoord([self.x,self.y])
		for i in range(0,len(self.shoppingList)):
			targetInd = store.getIndexFromCoord(self.shoppingList[i])   
			thisDist = store.staticGraph.shortest_paths_dijkstra(source=startInd,target=targetInd,weights=None)[0][0]
			if (thisDist < shortestDist):
				shortestDist = thisDist
				shortInd = i
		if shortInd is None:
			raise ValueError("No unblocked paths available for the customer!")
		else:
			self.shoppingList.insert(0,self.shoppingList.pop(shortInd))
			
	## helper function for checking if position on shopping list
	def itemFound(self):
		if not len(self.shoppingList):
			raise ValueError("list of targets empty!")

		itemPos=self.shoppingList[0]
		if self.x==itemPos[0] and self.y==itemPos[1]:
			return True
		return False 
	
			
	def spreadViralPlumes(self,store):
		sample = np.random.random()
		if (sample < self.probSpreadPlume) and not store.useDiffusion:
			store.plumes[self.x,self.y] = PLUMELIFETIME
		elif store.useDiffusion:
			##check if cough or just constant emission
			if (sample < self.probSpreadPlume):
				store.plumes[self.x,self.y] += PLUMECONCINC
				print "Customer coughed at ({},{})".format(self.x,self.y)
			else:
				store.plumes[self.x,self.y] += PLUMECONCCONT # according to 1 min of emission is same as 6 coughs

	## randomly add 1...N coordinates to the shoppping list 
	def initShoppingList(self,store,maxN):
		targetsDrawn = np.random.randint(maxN)+1
		while len(self.shoppingList)<targetsDrawn:
			tx = np.random.randint(store.Lx)
			ty = np.random.randint(1,store.Ly) #  running from 1 to avoid customers visiting the row of cashiers before leaving
			while store.blocked[tx,ty] or [tx,ty] in self.shoppingList or ([tx,ty] in store.exit) or (store.entrance[0]==tx and store.entrance[1]==ty) or (tx<1 or ty<1) or (tx<3 and ty<3) or not (store.blockedShelves[tx,ty-1] or (ty+1< store.Ly and store.blockedShelves[tx,ty+1]) or store.blockedShelves[tx-1,ty] or (tx+1< store.Lx and store.blockedShelves[tx+1,ty])): ##last ensures that customers avoid exits before leaving
				tx = np.random.randint(store.Lx)
				ty = np.random.randint(store.Ly)
			self.addTarget([tx,ty])
		self.initItemsList = len(self.shoppingList)
		self.cashierWaitingTime = CASHIERTIMEPERITEM*targetsDrawn
		return targetsDrawn

	## when customer leaves the store, return some statistics
	def getFinalStats(self):
		return self.infected,self.initItemsList,self.timeInStore,self.exposure, self.exposureTime, self.exposureTimeThres

	# v0.000 takes a totally random step
	def takeRandomStep(self, store):
		
		direction = np.random.permutation(len(DIRECTIONS))
		for i in range(len(direction)):
			step = DIRECTIONS[direction[i]]
			tmpPos = np.array([self.x, self.y], dtype=int)+step
			if tmpPos[0]<0 or tmpPos[0]>=store.Lx or tmpPos[1]<0 or tmpPos[1]>=store.Ly:
				continue
			elif store.blocked[tmpPos[0],tmpPos[1]]==1:
				continue
			else:
				store.blocked[self.x,self.y] = 0
				self.x = tmpPos[0]
				self.y = tmpPos[1]
				store.blocked[self.x,self.y] = 1
				break
		return self.x, self.y

	
	def atExit(self, store):
		for s in store.exit:
			if self.x==s[0] and self.y==s[1]:
				return 1
		return 0
		


	
class SmartCustomer(Customer):


	def takeStep(self, store):
		self.timeInStore +=1

		## check exposure and plume spreading here to make sure that this is done on every step
		if store.plumes[self.x,self.y] and not store.useDiffusion:
			self.exposure+=1
		elif store.plumes[self.x,self.y] and store.useDiffusion:
			self.exposure+=store.plumes[self.x,self.y]*store.dt
			if not self.infected:
				store.storeWideExposure+=store.plumes[self.x,self.y]*store.dt
			if store.plumes[self.x,self.y]>0: 
				self.exposureTime+=1
				if store.plumes[self.x,self.y]>EXPOSURELIMIT: 
					self.exposureTimeThres+=1
		if (self.infected):
			self.spreadViralPlumes(store)

		if self.waitingTime:
			self.waitingTime -=1 
			return self.x, self.y


		if not len(self.shoppingList):
			# head for exit
			if not self.atExit(store):
				self.shoppingList.append(store.getExit())
				self.headingForExit = 1
			##shopping list done, exit found
			elif self.atExit(store) and (self.cashierWaitingTime > 0):
				self.cashierWaitingTime -= 1
				return self.x, self.y
			else:
				store.blocked[self.x,self.y] = 0 #make sure exit is not blocked once customer leaves the simulation
				return -1,-1
		
		if self.itemFound():
			itemPos = self.shoppingList.pop(0)
			self.waitingTime = np.random.randint(MINWAITINGTIME,MAXWAITINGTIME)
			return itemPos
			
		## Get the path to the next target coordinate
		if self.path is None or not len(self.path):
			self.updateFirstTarget(store)	
			startInd = store.getIndexFromCoord([self.x,self.y])
			targetInd = store.getIndexFromCoord(self.shoppingList[0])    

			self.path = store.staticGraph.get_shortest_paths(startInd,to=targetInd)[0]
			self.path.pop(0)

			## sanity check, apparently target already found
			if not len(self.path):
				itemPos = self.shoppingList.pop(0)
				return itemPos

		## 	
		if not len(self.path):
			print self.x, self.y, self.shoppingList, self.headingForExit
 
		step = store.getCoordFromIndex(self.path[0])
		## check that the step is possible, i.e. no other customer blocking
		if not store.blocked[step[0], step[1]]:
			store.blocked[self.x,self.y] = 0
			self.x=step[0]
			self.y=step[1]
			store.blocked[self.x,self.y] = 1
			self.path.pop(0)
			if not len(self.path):
				self.path=None
		## sanity check for simulations with small, 'office-like' environment 
		elif store.Lx*store.Ly < 101 and self.timeInStore%180 == 0:
			store.createStaticGraph()
			self.path=None
		## if step blocked by customer, maybe take one random step and compute new route to the next item (during next step)
		elif (not self.headingForExit and np.random.rand()<BLOCKRANDOMSTEP) or np.random.rand()<BLOCKRANDOMSTEP*1e-2: ## no point in random step if queing for exit
			self.takeRandomStep(store)
			self.path=None


		return self.x, self.y
	


