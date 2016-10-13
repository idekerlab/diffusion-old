#!/usr/bin/env python

import networkx as nx
import re, math, os, sys, operator, random
from scipy.sparse import coo_matrix
# from scipy.linalg import expm
import imp
#import kernel_scipy
from optparse import OptionParser
import csv
import numpy as np
#schema=imp.load_source('ndexSchema','ndexSchema.py')
import ndexSchema as schema
#util=imp.load_source('ndexUtil','ndexUtil.py')
import ndexUtil as util
#nc=imp.load_source('ndexClient','ndexClient.py')
import ndex.client as nc
#kernel=imp.load_source('kernel_scipy','kernel_scipy.py')
import kernel_scipy as kernel
import operator

class TableEntityMapper:
    """This class is a container for a mapping table between two namespaces.
    Expects a filename for a file that is two tab delimited columns, source\ttarget
    Retutns a list (!) of ids"""

    def __init__(self, filename, missingStrList=[]):
        self.mappingTable = {}
        self.readFileToMappingTable(filename)
        self.missingStrList = missingStrList

    def readFileToMappingTable(self, filename):
        with open(filename, "r") as ins:
            for line in ins:
                arr = line.split('\t')
                if arr[0] not in self.mappingTable:
                    self.mappingTable[arr[0]] = [arr[1].rstrip()]
                else:
                    print("multiple mappings of %s detected" % arr[0])
                    self.mappingTable[arr[0]].append(arr[1])


    def map(self, entity):
        try:
            return self.mappingTable[entity]
        except:
            return self.missingStrList


def mapByPrefix(i, MGImapper=None, RGDmapper=None, HGNCmapper=None, targetPrefix='hgnc:'):
    st = i.lower()
    if st.startswith(targetPrefix):
        hgout = [i[len(targetPrefix):len(i)]]
    elif st.startswith('hgnc:') and HGNCmapper is not None:
        hgout = HGNCmapper.map(i[5:len(i)])
    elif st.startswith('mgi:') and MGImapper is not None:
        hgout = MGImapper.map(i[4:len(i)])
    elif st.startswith('rgd:') and RGDmapper is not None:
        hgout = RGDmapper.map(i[4:len(i)])
    else:
        hgout = [i]
    return hgout


def filterByPrefix(l, targetPrefix='hgnc:'):
    out = []
    for i in l:
        st = i.lower()
        if not st.startswith('bel:'):
            if st.startswith(targetPrefix):
                out.append(i[len(targetPrefix):len(i)])
            else:
                out.append(i)
    return out


class NdexToGeneSif(util.NetworkWrapper):
    def __init__(self, ndexNetwork, prefix='hgnc:', MGImapper=None, RGDmapper=None, HGNCmapper=None):
        self.prefix = prefix
        self.network = ndexNetwork
        self.supportToEdgeMap = {}
        self.citationToSupportMap = {}
        self.nodeLabelMap = {}
        self.termLabelMap = {}
        #new vars
        self.basetermMap = {}
        self.missingBaseterms = []
        #fix this guy for non-bel:
        self.nodeBaseterms = {}
        self.MGImapper = MGImapper
        self.RGDmapper = RGDmapper
        self.HGNCmapper = HGNCmapper

        for nodeId, node in ndexNetwork['nodes'].iteritems():
            self.nodeLabelMap[int(nodeId)] = self.getNodeLabel(node)
            self.nodeBaseterms[int(nodeId)] = self.getNodeBaseterms(node)

        for edge in ndexNetwork['edges'].values():
            for supportId in edge['supportIds']:
                supports = ndexNetwork['supports']
                support = supports[str(supportId)]
                if supportId in self.supportToEdgeMap:
                    edgeList = self.supportToEdgeMap[supportId]
                else:
                    edgeList = []
                edgeList.append(edge)
                self.supportToEdgeMap[supportId] = edgeList

        for supportId in self.supportToEdgeMap.keys():
            support = ndexNetwork['supports'][str(supportId)]
            citationId = support['citationId']
            if citationId in self.citationToSupportMap:
                supportIdList = self.citationToSupportMap[citationId]
            else:
                supportIdList = []
            supportIdList.append(supportId)
            self.citationToSupportMap[citationId] = supportIdList

    #this is new
    def edgeExpandToSIF(self, edge):
        """Takes an edge and simplifies to gene symbol and interaction type.
       All baseterms associated with subject inherit all interactions with all objects"""
        subjectBasetermList = []
        objectBasetermList = []
        predicateLabel = "missing"
        subjectId = edge['subjectId']
        objectId = edge['objectId']
        if subjectId in self.nodeBaseterms:
            subjectBasetermList=self.nodeBaseterms[subjectId]
        if objectId in self.nodeBaseterms:
            objectBasetermList=self.nodeBaseterms[objectId]

        predicateId = edge['predicateId']
        predicateLabel = util.stripPrefixes(self.getTermLabel(predicateId), targetPrefix=self.prefix)
        sifTriples = []

        subjMapped = []
        for subj in subjectBasetermList:
            subjMapped.extend(mapByPrefix(subj, self.MGImapper, self.RGDmapper, self.HGNCmapper, targetPrefix=self.prefix))
        subjMapped = filterByPrefix(subjMapped)

        objMapped = []
        for obj in objectBasetermList:
            objMapped.extend(mapByPrefix(obj, self.MGImapper, self.RGDmapper, self.HGNCmapper, targetPrefix=self.prefix))
        objMapped = filterByPrefix(objMapped)

        for subj in subjMapped:
            for obj in objMapped:
                sifTriples.append([subj, predicateLabel, obj])

        return sifTriples


    def getNodeBaseterms(self, node):
        if 'name' in node and node['name']:
            return [node['name']]

        elif 'represents' in node:
            return self.basetermMap[node['represents']]

        else:
            return ["node %s"] % (node['id'])


    def writeSIF(self, fileName=None, append=False):
        if fileName and not append:
            output = open(fileName, 'w')
        elif fileName and append:
            output = open(fileName, 'a')
        else:
            output = sys.stdout

        for edge in self.network['edges'].values():
            sifs = self.edgeExpandToSIF(edge)
            for s in sifs:
                if len(s) == 3:
                    output.write('\t'.join(s) + '\n')

        if fileName:
            output.close

    def edgeList(self):
        edge_list=list()

        for edge in self.network['edges'].values():
            sifs = self.edgeExpandToSIF(edge)
            for s in sifs:
                edge_list.append((s[0],s[2]))
        return edge_list

    def checkType(self,term):
        if 'name' in term.keys():
            return 'baseterm'
        elif 'functionTermId' in term.keys():
            return 'functionterm'
        elif 'edgeId' in term.keys():
            return 'reifiededgeterm'
        else:
            return 'unknown'

    def getTermLabel(self, termId):
        if termId in self.termLabelMap:
            return self.termLabelMap[termId]
        else:
            label = "error"
            term = self.getTermById(termId)
            #type = term['type'].lower()
            type=self.checkType(term)
            if type == "baseterm":
                name = term['name']
                if 'namespaceId' in term and term['namespaceId']:
                    namespaceId = term['namespaceId']
                    namespace = self.getNamespaceById(namespaceId)

                    if namespace:
                        if namespace['prefix']:
                            label = "%s:%s" % (namespace['prefix'], name)
                        elif namespace['uri']:
                            label = "%s%s" % (namespace['uri'], name)
                        else:
                            label = "%s%s" % ('BEL:', name)
                    else:
                        label = "%s%s" % ('BEL:', name)
                else:
                    label = "%s%s" % ('BEL:', name)
                self.basetermMap[termId] = [label]

            elif type == "functionterm":
                functionTermId = term['functionTermId']
                functionLabel = self.getTermLabel(functionTermId)
                functionLabel = util.getFunctionAbbreviation(functionLabel)
                parameterLabels = []
                if (termId not in self.basetermMap):
                    self.basetermMap[termId] = []
                for parameterId in term['parameterIds']:
                    parameterLabel = self.getTermLabel(parameterId)
                    parameterLabels.append(parameterLabel)
                    if parameterId in self.basetermMap:
                        self.basetermMap[termId].extend(self.basetermMap[parameterId])
                    else:
                        self.missingBaseterms.append(parameterId)
                label = "%s(%s)" % (functionLabel, ",".join(parameterLabels))

            elif type == "reifiededgeterm":
                edgeId = term['edgeId']
                edges = self.network['edges']
                self.basetermMap[termId] = []
                if edgeId in edges:
                    reifiedEdge = edges[edgeId]
                    label = "(%s)" % (self.getEdgeLabel(reifiedEdge))
                else:
                    label = "(reifiedEdge: %s)" % (edgeId)

            else:
                label = "term: %s" % (termId)
                self.basetermMap[termId] = []

            self.termLabelMap[termId] = label
            return label


def queryVector(query, labels):
    out = {}
    weight_sum = 0.0
    for i in labels:
        if i in query:
            out[i] = 1.0
            weight_sum += 1.0
        else:
            out[i] = 0.0

    for i in labels:
        out[i] = out[i] / weight_sum

    return out


def weightedQueryVector(query, labels, weights):
    out = {}
    weight_sum = 0
    for i in labels:
        if i in query:
            try:
                out[i] = weights[i]
                weight_sum += weights[i]
            except KeyError:
                out[i] = 1
                weight_sum += 1
        else:
            out[i] = 0

    for i in labels:
        out[i] = np.float32(out[i] / weight_sum)

    return out


def readNodeWeights(filename):
    weights = {}
    for l in open(filename, 'r'):
        line=l.rstrip().split('\t')
        if len(line)==2:
            weights[line[0]] = float(line[1])
        else:
            print "problem line:"+'TAB'.join(line)
    return weights


def readSif(filename):
    triples = []
    for line in csv.reader(open(filename, 'rb'), delimiter='\t'):
        triples.append(line)
    return triples


def filterSif(triples, scores, desired_nodes=30):
    scores_made_the_cut = sorted(scores.items(), key=operator.itemgetter(1), reverse=True)[0:(desired_nodes)]
    nodes_made_the_cut = [i[0] for i in scores_made_the_cut]
    edges_made_the_cut = []
    for tr in triples:
        if (tr[0] in nodes_made_the_cut) and (tr[2] in nodes_made_the_cut):
            edges_made_the_cut.append(tr)
    return edges_made_the_cut
