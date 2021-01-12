import numpy as np
import igraph as ig
import Image as im
from Params import *

np.set_printoptions(threshold=100000)



class Store:

	def __init__(self, Lx, Ly, dx):
		self.Lx = Lx
		self.Ly = Ly
		self.dt = 1.0
		self.dx =dx
		self.dy =dx 
		self.blocked = None	
		self.blockedShelves = None	
		self.invLxLy = 1./(Lx*Ly)
		self.staticGraph = ig.Graph()
		self.entrance=None
		self.useDiffusion=None
		self.exit=[]
		self.exitActive=np.zeros(NEXITS)
		self.plumes = np.zeros((self.Lx,self.Ly)) ##This contains natural numbers, 0 means no viral plume is present at the grid point,
		self.plumesNew = np.zeros((self.Lx,self.Ly)) ##This contains natural numbers, 0 means no viral plume is present at the grid point, while a value > 0 presents how many more temporal turns the plume will be present at a grid point
		self.plumesIntegrated = np.zeros((self.Lx,self.Ly)) ##This contains natural numbers, 0 means no viral plume is present at the grid point, while a value > 0 presents how many more temporal turns the plume will be present at a grid point
		self.diffusionCoeff = np.ones((self.Lx,self.Ly))*DIFFCOEFF
		self.ACSinkCoeff = np.ones((self.Lx,self.Ly))*ACSINKCOEFF ##Coefficient for the sink term of the form: -k*c.
		self.storeWideExposure = 0

	def initializeExposureDuringTimeStep(self):	
		self.storeWideExposure = 0

	## randomly place shelves into the store grid until area fraction af covered  
	def initializeShelves(self,af, maxSize):
		self.blocked = np.zeros((self.Lx,self.Ly))	
		runningAf = 0
		
		while runningAf < af:
			shelfSize=np.random.randint(2,maxSize,2)
			shelfPosx = np.random.randint(1,self.Lx-shelfSize[0]-1) 
			shelfPosy = np.random.randint(1,self.Ly-shelfSize[1]-1) 
			while self.blocked[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]].sum():
				shelfPosx = np.random.randint(1,self.Lx-shelfSize[0]-1) 
				shelfPosy = np.random.randint(1,self.Ly-shelfSize[1]-1) 
			self.blocked[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]] = 1
			self.diffusionCoeff[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]] = 0
			runningAf += shelfSize[0]*shelfSize[1] * self.invLxLy

		return runningAf

	## initialize the store layout from a separate file	
	def loadImageAsGeometry(self, filename):
		img = im.open(filename)
		#~ img.show()
		px = img.load()
		#~ for some reason as a result of the loading routine, the subsequent px has 0=black and 255=white, i.e. forces a pixel map from 0-255
		self.blocked = np.zeros((self.Lx,self.Ly))	
		self.free = np.zeros((self.Lx,self.Ly))	
		for i in range(0,self.Lx):
			for j in range(0,self.Ly):
				if not px[i,j]:
					self.blocked[i,j] = 1
					self.diffusionCoeff[i,j] = 0
		self.blockedShelves=np.copy(self.blocked)

					

	def save_image(self, npdata, outfilename ):
		img = im.fromarray(np.asarray(np.clip(npdata,0,1), dtype="int"), "L")
		img.save(outfilename)
	
	def updateDiffusion(self):
		self.plumesNew = np.copy(self.plumes)
		self.plumesNew[1:,:]  +=  self.diffusionCoeff[1:,:]  * self.diffusionCoeff[:-1,:] * self.dt * (self.plumes[:-1,:] - self.plumes[1:,:])/(self.dx**2) /DIFFCOEFF 
		self.plumesNew[:-1,:]  += self.diffusionCoeff[:-1,:] * self.diffusionCoeff[1:,:]  * self.dt * (self.plumes[1:,:] - self.plumes[:-1,:])/(self.dx**2) /DIFFCOEFF
		self.plumesNew[:,1:]  +=  self.diffusionCoeff[:,1:]  * self.diffusionCoeff[:,:-1] * self.dt * (self.plumes[:,:-1] - self.plumes[:,1:])/(self.dy**2) /DIFFCOEFF
		self.plumesNew[:,:-1]  += self.diffusionCoeff[:,:-1] * self.diffusionCoeff[:,1:]  * self.dt * (self.plumes[:,1:] - self.plumes[:,:-1])/(self.dy**2) /DIFFCOEFF
		
		self.plumesNew[:,:] -= self.ACSinkCoeff[:,:] * self.dt * self.plumes[:,:] #Finally, apply the sink term where applicable
	
		self.plumes = self.plumesNew
		self.plumes[self.plumes<PLUMEMIN] = 0

		## compute the plume concentration through siumulation
		self.plumesIntegrated += self.plumes
		return
	
	## add discrete plumes  NOTE: not in use with diffusion
	def addPlume(self, plumeDuration):
		plumePosx = np.random.randint(1,self.Lx-1)
		plumePosy = np.random.randint(1,self.Ly-1)
		while self.blocked[plumePosx,plumePosy] or self.plumes[plumePosx, plumePosy]:
			plumePosx = np.random.randint(1,self.Lx-1)
			plumePosy = np.random.randint(1,self.Lx-1)
		self.plumes[plumePosx,plumePosy] = plumeDuration
		return  np.array([plumePosx, plumePosy])	

	## if runnning simulation with one customer, add static field of plumes to the store	
	def initStaticPlumeField(self, nPlumes):
		for i in range(nPlumes):
			self.addPlume(1)
		return
		
	
	## initialize the store layout with semi-regular grid:  place N shelves of size I length and J width with distance D
	def initializeShelvesRegular(self, N, I=2, J=1, D=2 ): ## dimensions in meters
		self.blocked = np.zeros((self.Lx,self.Ly))	
		placed = 0
		tries = 0
		II = int(I/self.dx)
		JJ = int(J/self.dx)
		DD = int(D/self.dx)

		while placed<N:
			## randomly in x or y-direction
			if np.random.rand()<0.5:
				shelfSize=np.array([II,JJ])
				axis = 1
			else:
				shelfSize=np.array([JJ,II])
				axis =0
			shelfPosx = np.random.randint(1,self.Lx-shelfSize[0]-1)
			shelfPosy = np.random.randint(1,self.Ly-shelfSize[1]-1)
		
			while self.blocked[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]].sum():
				shelfPosx = np.random.randint(1,self.Lx-shelfSize[0]-1) 
				shelfPosy = np.random.randint(1,self.Ly-shelfSize[1]-1) 
				tries+=1
				if tries>1e4:
					self.blockedShelves=np.copy(self.blocked)
					return placed
			self.blocked[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]] = 1
			placed+=1			

			direction = np.random.choice([-1,1])
			while placed<N:
				if axis:
					shelfPosy+=direction*(DD+JJ)
				else:
					shelfPosx+=direction*(DD+II)
				if shelfPosx<0 or shelfPosx>=self.Lx-shelfSize[0]-1 or shelfPosy<0 or shelfPosy>=self.Ly-shelfSize[1]-1:
					break 
				if not self.blocked[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]].sum():
					self.blocked[shelfPosx:shelfPosx+shelfSize[0],shelfPosy:shelfPosy+shelfSize[1]] = 1
					placed+=1
				else: 
					break
		self.blockedShelves=np.copy(self.blocked)
		return placed	
				



	## compute the navigation graph for the store layout	
	def createStaticGraph(self):
		totNodes = self.Lx*self.Ly
		blockedNodeList = np.asarray(self.blocked).reshape(-1)
		self.staticGraph.add_vertices(totNodes)
		#~ Connect the nearest neightbors along the x-axis
		for i in range(0,totNodes):
			try:
				if (blockedNodeList[i]==0) and (blockedNodeList[i+1] == 0) and (self.staticGraph.are_connected(i,i+1) == False) and ((i+1) % self.Lx != 0):
					self.staticGraph.add_edges([(i,i+1)])
			except:
				continue
				
		#~ Connect the nearest neighbors along the y-axis
		for i in range(0,totNodes):
			try:
				if (blockedNodeList[i]==0) and (blockedNodeList[i+self.Lx] == 0) and (self.staticGraph.are_connected(i,i+self.Lx) == False):
					self.staticGraph.add_edges([(i,i+self.Lx)])
			except:
				continue

				
	## helper functions to get from index to coordinate and vice versa for the route planning algorithm
	def getCoordFromIndex(self, idx):
		return [int(np.floor(idx/self.Lx)), idx%self.Lx ]

	def getIndexFromCoord(self,coord):	
		return coord[1]+self.Lx*coord[0]


	## randomly place doors to the store 
	def initializeDoors(self):
		## by default, exit and entrance are in corners
		entranceInd = self.getIndexFromCoord([ENTRANCEPOS,0])
		self.entrance= self.getCoordFromIndex(entranceInd)


		i = EXITPOS
		while len(self.exit)<NEXITS:
			exitInd = self.getIndexFromCoord([self.Lx-i,0])

			## ensure that there is a path between entrance and exit for obvious reasons
			checkPossiblePath = self.staticGraph.shortest_paths_dijkstra(source=entranceInd,target=exitInd,weights=None)[0][0]
			while np.isinf(checkPossiblePath) and self.Lx-i>0:
				i+=1	
				exitInd = self.getIndexFromCoord([self.Lx-i,0])
				checkPossiblePath = self.staticGraph.shortest_paths_dijkstra(source=entranceInd,target=exitInd,weights=None)[0][0]
			
			self.exit.append(self.getCoordFromIndex(exitInd))
			i+=CASHIERD
		## sanity check
		print "exits at : ", self.exit 
		if np.isinf(checkPossiblePath):
			raise ValueError("no path between store entrance and exit!")
		self.exit = self.exit[::-1]
		return	

	## method to give a cashoier for the leaving customer; the exit with the least customers is given
	def getExit(self):
		exitInd = np.argmin(self.exitActive)
		self.exitActive[exitInd]+=1
		return self.exit[exitInd]

	def updateQueue(self, exitPos):
		self.exitActive[self.exit.index(exitPos)]-=1
		return


