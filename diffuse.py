__author__ = 'decarlin'

import sys
import logging
import operator
import json
import uuid
import csv

import ndex.client as nc

import diffusiond.src.kernel_scipy as kernel
import diffusiond.src.automateHeatKernel as ahk

from diffusiond.src.multiplyEstablishedKernel import EstablishedKernel

def diffuse(cx, identifier_set):
    """Diffuses a network represented as cx against an identifier_set"""
    logging.info('Diffusion is now converting the CX json to the SIF format')
    ndex_network = cxToNDEx(cx)
    sif_network = ahk.Ndex.ToGeneSif(ndex_network)
    kernel_id = create_kernel(sif_network)
    ranked_entities = diffuse_against_kernel(kernel_id, identifier_set)
    logging.info('Diffusion completed, now returning the ranked entities as json to the caller')
    return json.dumps(ranked_entities)

def create_kernel(sif_network):
    """Generate a kernel for a Sif network and write it to disk, returning a uuid."""
    logging.info('Diffusion is now creating a kernel for the sif network')
    edges = sif_network.edgeList()
    diffusion_kernel = kernel.SciPyKernel(edges)
    kernel_id=uuid.uuid1()
    logging.info('Diffusion is now writing the kernel with kernel id: {0} to disk'.format(str(kernel_id)))
    diffusion_kernel.writeKernel('kernels/{0}'.format(str(kernel_id)))
    return kernel_id

def diffuse_against_kernel(kernel_id, identifier_set):
    """Diffuse an identifier_set against the generated kernel"""
    logging.info('Diffusion is now loading the kernel with kernel id: {0} for use'.format(str(kernel_id)))
    ker = EstablishedKernel('kernels/{0}'.format(kernel_id))
    logging.info('Diffusion is now diffusing against the kernel with the identifier_set')
    queryVector = ahk.queryVector(identifier_set, ker.labels)
    diffused_network=ker.diffuse(queryVector)
    return sorted(diffused.items(), key=operator.itemgetter(1), reverse=True)

def cyToNDEx(cx):
  """Converts cx json into NDEx json"""
  #TODO IMPL
  return cx

def main():
  args = sys.argv
  num_args = len(args)
  logging.info('Diffusion was called with {0} arguments'.format(num_args))
  if num_args != 2:
    error_message = 'Diffusion expected 2 arguments, recieved: {0}\n'.format(num_args)
    log.error(error_message)
    sys.stderr.write(error_message)
  else:
    log.info('Evoking diffuse...')
    ranked_entities = diffuse(args[1], args[2])
    log.info('Writing reply...')
    print(ranked_entites)


