#!/usr/bin/env python
from __future__ import division
import networkx as nx
import csv
import matplotlib.pylab as plt
import random as rand
import cmath
from operator import itemgetter
from numarray import *
 
#construct the Infocom graph from graph.dat file
def buildG_info( G ):
	reader = csv.reader(file("/home/kazem/Desktop/traces/mobility/graph_inf06.dat"))
	
	for line in reader:
		G.add_edge(line[0],line[1],weight=float(line[2]))
	
#draw the contact graph
def draw_Graph( G ):
	nx.draw_spring(G,
			node_color=[float(G.degree(v)) for v in G],
			node_size=40,
			with_labels=False,
			cmap=plt.cm.Reds,
			)
	plt.show()

#this generates a random external list with size of num members	
def gen_Ext( G, num ):	
	n_ = []	#nodes list
	for i in G.nodes():
		n_.append(i)
	e_ = []	#external list
	e_ = rand.sample(n_, num)
	return e_
#this method deletes all edges among external nodes in G				
def rem_Ext( G, ext_ ): 	
	for e_ in G.edges():
		pairs_ = []
		for i in e_:
			pairs_.append(i)
		if(ext_.count(pairs_[0]) == 1 and  ext_.count(pairs_[1]) == 1):
			G.remove_edge(pairs_[0],pairs_[1])
		#print pairs_[0]
		#print pairs_[1]

#generate the sorted list of edges among external node and store it in a dictionary
def gen_base( G_ ):
	base_ = []
	for e_ in G_.edges():
		pairs_ = []	
		for i in e_:
			pairs_.append(i)
		#print "%s %s %f" % (pairs_[0], pairs_[1], G_[pairs_[0]][pairs_[1]]['weight'])
		if(int(pairs_[0]) <= int(pairs_[1])):
			e_ = ( pairs_[0], pairs_[1], G_[pairs_[0]][pairs_[1]]['weight'] )
			obj_ = ()
			for tup_ in base_:
				if (tup_[2] >= G_[pairs_[0]][pairs_[1]]['weight'] ):
					obj_ = (tup_[0], tup_[1], tup_[2])
					break
			if(obj_ != ()):
				i = base_.index(obj_)
				base_.insert(i, e_)	
			else:
				base_.append( e_ )
		else:
			e_ = (pairs_[1],pairs_[0], G_[pairs_[0]][pairs_[1]]['weight'])
			obj_ = ()
			for tup_ in base_:
				if (tup_[2] >= G_[pairs_[0]][pairs_[1]]['weight'] ):
					obj_ = (tup_[0], tup_[1], tup_[2])
					break
			if(obj_ != ()):
				i = base_.index(obj_)
				base_.insert(i, e_)	
			else:
				base_.append( e_ )
	return base_

#this method computes predictions based on common number of neighbors among external nodes		
def predCN( G, G_, ext_):
	neigh_ = {}	# a dictionary for keeping neighbor set of all external nodes
	predCN_ = []	# a list for predictions based on number of common neighbors
	for u in ext_:
		#let's find neighbor set of u
		cnt = 0
		inters = 0
		if u not in neigh_:
			for i in G.neighbors(u):
				if(cnt == 0):
					neigh_[u] = set(i)
					cnt = 1
				else:
					neigh_[u].add(i)
			
		for v in ext_:
			cnt = 0
			if v not in neigh_:
				for i in G.neighbors(v):
					if(cnt == 0):
						neigh_[v] = set(i)
						cnt = 1
					else:
						neigh_[v].add(i)
			if(u < v):
				if(G_.has_edge(u, v)):
					#we add a score if the corresponding edge is in G_
					#print "(%s, %s)" % (u,v)
					#print neigh_[u]
					#print neigh_[v]
					if (u in neigh_ and v in neigh_):
						inters = len(neigh_[u].intersection(neigh_[v]))
						e_ = (u, v, inters)
					else:
						inters = 0
						e_ = (u, v, inters)
					# we have to find the place we should add this pair based on theor score
					obj_ = ()
					for tup_ in predCN_:
						if (tup_[2] <= inters):
							obj_ = (tup_[0], tup_[1], tup_[2])
							break
					if(obj_ != ()):
						i = predCN_.index(obj_)
						predCN_.insert(i, e_)	
					else:
						predCN_.append( e_ )
					#print u,v,len(neigh_[u].intersection(neigh_[v]))
	return predCN_
									
#this method computes predictions based on SPF among external nodes
def predSPF( G, G_, ext_):			
	length = nx.shortest_path_length(G,weighted=True)
	#print length
	predSPF_ = []	# a list for predictions based on SPF
	for u in ext_:
		for v in ext_:
			if(u < v):
				if(G_.has_edge(u, v)):
					if (u in length and v in length[u]):
						len_ = length[u][v]
						e_ = (u, v, len_)
					else:
						len_ = 1000000.0
						e_ = (u, v, len_)
						
					obj_ = ()
					for tup_ in predSPF_:
						if (tup_[2] >= len_):
							obj_ = (tup_[0], tup_[1], tup_[2])
							break
					if(obj_ != ()):
						i = predSPF_.index(obj_)
						predSPF_.insert(i, e_)	
					else:
						predSPF_.append( e_ )
					#print u,v,
	return predSPF_
	
def evalPreds( G_, predL, extRank_, msg):	
	results_ = zeros([1,11], Float)
	inx = 0 
	#for limit in range(0, G_.number_of_edges(), G_.number_of_edges()//10):
	limit = 0.0
	while limit <= G_.number_of_edges():
		#print limit
		cnt1 = 0;
		cnt2 = 0;
		match = 0
		for guess_ in predL:
			cnt1 = cnt1 + 1
			cnt2 = 0
			if(cnt1 <= limit):
				#let's search this item in external list
				for ex_ in extRank_:
					cnt2 = cnt2 + 1	
					if(cnt2 <= limit):
						if((guess_[0] == ex_[0] and guess_[1] == ex_[1]) or (guess_[0] == ex_[1] and guess_[1] == ex_[0])):
							match = match + 1
							#print guess_
					else:
						break
			else:
				break
		if(limit > 0):
			#results_.append(match/limit)
			results_[0][inx] = match/limit
			inx += 1
		else:
			#results_.append(0.0)
			results_[0][inx] = 0.0
			inx += 1
		#if(G_.number_of_edges() > 0):
		#	rand_.append(limit/G_.number_of_edges())
		#print "%d matches from %d" % (match_,limit)
		
		limit += G_.number_of_edges()/10.0
		#print limit	
	return results_
	#print "%s predictor results:" % msg	
	#print results_
	
	#print "random results:"	
	#print rand_
	
#method: print the prediction results in the list
def printList( predList ):
	for e_ in predList:
		print e_
#this method simply picks a neighbor of node u in G based on their transition probabilities
def random_pick(graph, u):
	x = rand.uniform(0,1)
	cumulative_probability = 0.0
	for v in graph.neighbors(u):
		cumulative_probability += graph[u][v]['weight']
		if x < cumulative_probability: break
	return v
#a method which runs a random walk		
def unweighted_random_walk(graph, starting_point, ext_, T, visiting_no, alpha):
	
	'''
	starting_point: String that represents the starting point in the graph
	ending_point: String that represents the ending point in the graph
	graph: A NetworkX Graph object
	'''
	##Begin the random walk
	current_point = starting_point
	#visiting_no = 0
	steps = 0
	#Determine the hitting time to get to an arbitrary neighbor of the
	#starting point
	while steps <= T:
		steps += 1
		x = rand.uniform(0,1)
		if (x <= alpha):
			#return to the starting node with probability alpha
			current_point = starting_point
		else:
			#move to a random neighbor of current node with probability of 1-alpha
			#pick one of the edges out of the starting_node with equal probs
			next = random_pick(graph, current_point)
			#print "next node: %s" % next
			current_point = next
		#if ( current_point == ending_point ):
		if ( ext_.count(current_point) == 1 and current_point != starting_point ):
			visiting_no[ext_.index(starting_point)][ext_.index(current_point)] += 1
			#visiting_no += 1
	return visiting_no
		
#contact graph
spf_ = zeros([1,11], Float)
cn_ = zeros([1,11], Float)
rnd_ = zeros([1,11], Float)
prank_ = zeros([1,11], Float)
times = 0
for i in range(times):
	G = nx.Graph()
	G_ = nx.Graph()	#the graph induced by external nodes
	#build the contact graph from input file
	buildG_info( G )
	
	#let's generate a random external list with size of next_
	nExt_ = 0	#no of external nodes
	ext_ = []	#random external list
	ext_ = gen_Ext( G, nExt_ )
	#let's simulate the external nodes
	#ext_ = ['26', '33', '36', '41', '44', '45', '46', '55', '68', '80']
	G_ = G.subgraph(ext_)	#graph induced by external nodes
	
	#let's remove all edges among external nodes from the contact graph
	print "no of edges before removal: %d" % G.number_of_edges()
	rem_Ext(G, ext_)
	print "no of edges after removal: %d" % G.number_of_edges()
	
	#let's generate all edges among externals and sort them for posporessing
	extRank_ = []
	extRank_ = gen_base( G_ )
	
	#random predictor
	predRand_ = []	#let's shuffle this as a random predictor
	predRand_ = gen_base( G_ )
	rand.shuffle(predRand_)
	
	ext_.sort()
	draw_Graph( G )
	
	#print "non-sorted external nodes:"
	#for i in ext_:
		#print i	
	
	#generate predictions based on CNN
	predCN_ = predCN( G, G_, ext_)
	
	#generate predictions based on SPF
	predSPF_ = predSPF(G, G_, ext_)	
	
	#let's find total number of edges among external nodes
	print "no of edges in G_: %d" % G_.number_of_edges()
	
	#let's print external nodes edges with their weights
	#printList( extRank_ )	
	
	#let's print our CN predictor results
	#printList( predCN_ )	
	
	#let's print SPF predictor results
	#printList( predSPF_ )
		
	#let's evaluate CN predictor
	cn_ += evalPreds( G_, predCN_, extRank_, 'CN')
	
	#let's evaluate SPF predictor
	spf_ += evalPreds( G_, predSPF_, extRank_, 'SPF')
	
	rnd_ += evalPreds( G_, predRand_, extRank_, 'Random')
	
	T = 10000
	#Gt = nx.Graph()
	#Gt.add_edge(1,2,weight=1)
	#Gt.add_edge(1,3,weight=2)
	#Gt.add_edge(2,3,weight=3)
	#first let's generate the G's equivalent directed one
	#G.clear()
	#G = nx.complete_graph(5)
	#for u in G.nodes():
		#for v in G.neighbors(u):
			#G[u][v]['weight'] = 1
	#map the undirected graph G to its equivalent directed one H		
	H = G.to_directed()
	#remove edges among external nodes in H
	#ext_ = [0, 1, 4]
	rem_Ext(H, ext_)
	#print H.edges()
	#let's translate weights of G to probabilities in H
	for u in G.nodes():
		sum = 0.0
		for w in G.neighbors(u):
			sum += 1.0/G[u][w]['weight']
		for v in G.neighbors(u):
			H[u][v]['weight'] =  1.0/(G[u][v]['weight']*sum)
	
	visiting_no = zeros([nExt_,nExt_], Float)	##visiting matrix
	#print "a walk starting from %s in %d steps" % ( ext_[0], T )
	#for w in H.neighbors(ext_[0]):
	#	print "node: %s, p() = %f" % (w,H[ext_[0]][w]['weight'])
	alpha = 0.1
	for u in ext_:
		visits = unweighted_random_walk(H, u, ext_, T, visiting_no, alpha )
		#print visits
	#print "no of visists: %d, and p(.)=%f" % (visits[0][3],visits[0][3]/T)
	
	visits = visits/T
	#print visits
	predPRanks_ = []	# a list for predictions based on Page Rank
	for u in ext_:
		for v in ext_:
			if(u < v):
				if(G_.has_edge(u, v)):
					pr_ = visits[ext_.index(u)][ext_.index(v)] + visits[ext_.index(v)][ext_.index(u)]
					e_ = (u, v, pr_)
					obj_ = ()
					for tup_ in predPRanks_:
						if (tup_[2] <= pr_):
							obj_ = (tup_[0], tup_[1], tup_[2])
							break
					if(obj_ != ()):
						i = predPRanks_.index(obj_)
						predPRanks_.insert(i, e_)	
					else:
						predPRanks_.append( e_ )
					#print u,v,
	#let's print SPF predictor results
	#printList( predPRanks_ )
	
	#let's evaluate PageRank predictor
	prank_ += evalPreds( G_, predPRanks_, extRank_, 'PR')
if(times > 0):
	cn_ /= times 
	spf_ /= times
	prank_ /= times
	rnd_ /= times
	
	print "CN predictor results:"	
	print cn_
	
	print "SPF predictor results:"	
	print spf_
	
	print "PageRank predictor results:"	
	print prank_
	
	print "Random predictor results:"	
	print rnd_

G = nx.Graph()
#build the contact graph from input file
buildG_info( G )
draw_Graph( G )
length = nx.shortest_path_length(G,weighted=True)
#print length
u = '22'
#for v in G.nodes():
	#if (u in length and v in length[u]):
		#print length[u][v]
print length['35']
items=length['35'].items()
backitems=[ [v[1],v[0]] for v in items]
backitems.sort()
sortedlist=[ backitems[i][1] for i in range(0,len(backitems))]
print sortedlist		