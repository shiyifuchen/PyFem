# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 19:51:17 2014

@author: LUX
"""

from pyfem.util.dataStructures import Properties

def containsValue( db , val ):

  keys = db.keys()
  
  for key in keys:

    if type(db[key]) == dict:
      if containsValue(db[key],val):
        return True
    else:
      if db[key] == val:
        return True
  return False

def getType( a ):

  if a == 'true':
    return True
  elif a == 'false':
    return False
  else:
    try:
      return eval(a)
    except:
      return a

def storeValue( db , key , a ):

  if type(a) == list:
    tmp=[]
    for v in a:
      tmp.append(getType(v))
    db.store( key , tmp )
  else:
    db.store( key , getType(a) )

def readItem( l1 , db ):

  if '.' in l1[0]:
    l2 = l1[0].split('.',1)

    if db.has_key(l2[0]):
      if type(db[l2[0]]) == dict:
        child=db[l2[0]]
      else:
        child=Properties()
    else:
      child=Properties()

    l1[0]=l2[1]

    ln = readItem( l1 , child )

    db[l2[0]] = child

    return ln

  else:
    l2 = l1[1].split(';',1)

    if l2[0][0] == '[':
      l3 = l2[0][1:-1].split(',')
      storeValue( db , l1[0] , l3 )
    else:
      storeValue( db , l1[0] , l2[0] )
   
    return l2[1]

def readBlock( ln , db ):

  while True:
    #处理需要include其他文件的情况      
    if ln[0:7] == 'include':
      l1 = ln.split(';',1)
      deepFileParser( l1[0][8:-1] , db )
      ln = l1[1]
      continue

    l1 = ln.split('=',1)

    if len(l1) == 1:
      return ln

    if l1[0][0:2] == '};':
      return ln[2:]

    if l1[0][0:2] == '//':
      ln = l1[1].split(';',1)[1]
      continue

    if l1[1][0] == '{':
      child = Properties()
      ln = l1[1][1:]
     
      ln = readBlock( ln , child )

      db.store( l1[0] , child )

    else:     
      ln = readItem( l1 , db )

def fileParser( fileName ):

  db = Properties()

  f = open(fileName)
  
  f2 = ''
 
  for line in f:
    if not line.startswith('#'):
      f2 = f2+line
  #去除注释行
    
  ln = f2.replace('\n','').replace('\t','').replace(' ','').replace('\r','')

  readBlock( ln , db )

  return db

#处理include的文件
#此处未考虑有注释的情况，可修改
def deepFileParser( fileName , db ):

  ln = open(fileName).read().replace('\n','').replace('\t','').replace(' ','').replace('\r','')

  readBlock( ln , db )

  return db