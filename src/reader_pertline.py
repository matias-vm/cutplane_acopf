#!/usr/bin/python

import sys
import math
import cmath
import numpy as np
from random import randint

from myutils import *
import time
from log import *


class Bus:
  def __init__(self, count, nodeID, nodetype, Pd, Qd, Gs, Bs, Vbase, Vmax, Vmin, busline0):
    self.count = count
    self.nodeID = nodeID
    self.nodetype = nodetype
    self.Pd = Pd
    self.Qd = Qd
    self.genidsbycount = []
    self.frombranchids = {}
    self.tobranchids = {}
    self.Gs = Gs
    self.Bs = Bs
    self.Vbase = Vbase
    self.Vmax = Vmax
    self.Vmin = Vmin
    self.Pbalance = 0
    self.outdegree = self.indegree = self.degree = 0
    self.busline0 = busline0

    self.lat = -1
    self.lon = -1


  def getbusline0(self):
    return self.busline0
  def addgenerator(self, log, generatorcount, generator):
    self.genidsbycount.append(generatorcount)
    loud = 0 
    if loud:
      log.joint(" added generator # " + str(generatorcount) + " to bus ID " + str(self.nodeID))
      log.joint(" Pmax " + str(generator.Pmax) + " Pmin " + str(generator.Pmin) + "\n")

  def addfrombranch(self, log, id):
    quant = len(self.frombranchids)
    self.frombranchids[quant] = id
    self.outdegree += 1
    self.degree += 1
  def addtobranch(self, log, id):
    quant = len(self.tobranchids)
    self.tobranchids[quant] = id
    self.indegree += 1
    self.degree += 1


class branch:
    def __init__(self, log, count, f, id_f, t, id_t, r, x, bc, rateAmva, rateBmva, rateCmva, ratio, angle, maxangle, minangle, status, defaultlimit, branchline0):
       self.count = count
       self.f = f
       self.t = t
       self.id_f = id_f
       self.id_t = id_t
       self.r = r
       self.x = x
       self.bc = bc
       self.count = count
       self.branchline0 = branchline0
       self.rateAmva = rateAmva
       self.rateBmva = rateBmva
       self.limit = rateAmva
       self.constrainedflow = 1
       if self.limit == 0:
         # log.joint(" setting limit for branch " + str(count) + " " + str(defaultlimit) + "\n")
         self.limit = defaultlimit
         self.constrainedflow = 0
       if ratio == 0:
         ratio = 1
       self.ratio = ratio
       self.angle = angle
       self.angle_rad = math.pi*angle/180.0
       self.maxangle = maxangle
       self.maxangle_rad = math.pi*maxangle/180.0
       self.minangle = minangle
       self.minangle_rad = math.pi*minangle/180.0

       self.upperanglenone = 0
       if maxangle == 360 or maxangle == 0:
         self.upperanglenone = 1

       self.loweranglenone = 0         
       if minangle == -360 or minangle == 0:
         self.loweranglenone = 1

       self.invratio2 = invratio2 = 1/ratio**2
       self.multtf = multtf = 1/(ratio*cmath.exp(1j*self.angle_rad))
       self.multft = multft = 1/(ratio*cmath.exp(-1j*self.angle_rad))
       #print 'multtf', multtf
       self.status = status
       pi = math.pi
       self.z = z = r + x*1j
       self.y = y = 1/z
       self.Yff = Yff = (y + bc/2*1j)*invratio2
       self.Yft = Yft = -y*multft
       self.Ytf = Ytf = -y*multtf
       self.Ytt = Ytt = y + bc/2*1j
       self.Gff = Gff = (self.Yff).real
       self.Bff = Bff = (self.Yff).imag
       self.Gft = Gft = (self.Yft).real
       self.Bft = Bft = (self.Yft).imag
       self.Gtf = Gtf = (self.Ytf).real
       self.Btf = Btf = (self.Ytf).imag
       self.Gtt = Gtt = (self.Ytt).real
       self.Btt = Btt = (self.Ytt).imag
       loud = 0
#       if self.angle_rad != 0:
#         loud = 1
       if loud:
         log.joint("\nbr " + str(count) + " f " + str(f) + " t " + str(t) +"\n")
         log.joint("   idf " +  str(id_f) + " idt " + str(id_t) +"\n")
         log.joint("   r " + str(r) + " x " + str(x) + " bb " + str(bc) +"\n")
         log.joint("   ratio " + str(self.ratio) + " angle " + str(angle) + " angle_rad: " + str(self.angle_rad) + "\n")
         log.joint("   y " + str(y) + "\n")
         log.joint("       Yff " + str(Yff) + " , Yft " + str(Yft) + " , Ytf " + str(Ytf) +  " , Ytt " + str(Ytt) + "\n")

    def getbranchline0(self):
      return self.branchline0
    def show(self, log):
      log.joint(" < " + str(self.f) + " , " + str(self.t) + " > ")
      log.joint(" r " + str(self.r) + " x " + str(self.x) + " bc " + str(self.bc))
      log.joint("\n")
      log.joint(" ra " + str(self.ratio) + " ang " + str(self.angle) )
      log.joint("\n")
class gen:
    def __init__(self, count, nodeID, Pg, Qg, status, Pmax, Pmin, Qmax, Qmin, line0):
       self.count = count
       self.nodeID = nodeID
       self.Pg = Pg
       self.Qg = Qg
       self.status = status
       self.Pmax = Pmax
       self.Pmin = Pmin
       self.Qmax = Qmax
       self.Qmin = Qmin
       self.line0 = line0
       self.costlinenum = -1
    def addcost(self, log, costvector, linenum):
      self.costvector = costvector
      self.costdegree = len(costvector) - 1
      self.costlinenum = linenum
    def getline0(self):
      return self.line0


def readcase(log, all_data, casefilename):    
    log.joint("reading case file " + casefilename + "\n")

    t0 = time.time()

    try:
        f = open(casefilename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file " + casefilename)
        sys.exit("failure")

    readcode = readcase_thrulines(log, all_data, lines)
    all_data['casefilelines'] = lines

    t1 = time.time()

    log.joint("read time: " + str(t1 - t0) + "\n")



def readcase_thrulines(log, all_data, lines):
    log.joint("reading case file thru lines\n")

    pertline = 1
    ber      = 0
    normal   = 0
    normal0  = 1
    normal5  = 0

    if pertline:
      topfgname = 'top_flows_gens/' + all_data['casename'] + '_fg.txt' 

      log.joint("reading top flows and gens file " + topfgname + "\n")

      try:
        topfg    = open(topfgname,"r")
        topfg_lines    = topfg.readlines()
        topfg.close()
      except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

      line_top_branch   = topfg_lines[5].split()
      line_top_actgen   = topfg_lines[12].split()
      line_top_reactgen = topfg_lines[23].split()

      topbranch         = int(line_top_branch[1])
      topbranch_f       = int(line_top_branch[3])
      topbranch_t       = int(line_top_branch[5])
      topbranch_flow    = float(line_top_branch[5])
      topbranch_fdeg    = int(line_top_branch[9])
      topbranch_tdeg    = int(line_top_branch[11])

      log.joint(' the topbranch is ' + str(topbranch) + '\n')
      breakexit('c')

      if ((topbranch_flow < - 1e-8) and topbranch_fdeg == 1):
        log.joint(' perturbed instance should be infeasible\n')
        log.joint(' topbranch ' + str(topbranch) + ' flow ' + str(topbranch_flow) + ' f ' + str(topbranch_f) + ' degree ' + str(topbranch_fdeg) + '\n')
      elif ((topbranch_flow > 1e-8) and topbranch_tdeg == 1):
        log.joint(' perturbed instance should be infeasible\n')
        log.joint(' topbranch ' + str(topbranch) + ' flow ' + str(topbranch_flow) + ' f ' + str(topbranch_t) + ' degree ' + str(topbranch_tdeg) + '\n')
      
      
    newcasename   = '../data/' + all_data['casename'] + '_line5.m' 
    #newcasename   = '../data/' + all_data['casename'] + '_line.m' 
    #newcasename   = '../data/' + all_data['casename'] + '_b.m' 
    #newcasename   = '../data/' + all_data['casename'] + '_n_0_5.m' 
    #newcasename   = '../data/' + all_data['casename'] + '_n_5_5.m' 
    
    newcase       = open(newcasename,"w+")
    lookingforbus = 1
    linenum       = 2
    numlines      = len(lines)
    numisolated   = 0

    logpcasename  = 'perturbedcases_logs/' + all_data['casename'] + '_line5.log' #here we keep track of perturbations
    #logpcasename  = 'perturbedcases_logs/' + all_data['casename'] + '_line.log' #here we keep track of perturbations
    #logpcasename  = 'newcases_logs/' + all_data['casename'] + '_b.log' #here we keep track of perturbations
    #logpcasename  = 'newcases_logs/' + all_data['casename'] + '_n_0_5.log' #here we keep track of perturbations
    #logpcasename  = 'newcases_logs/' + all_data['casename'] + '_n_5_5.log' #here we keep track of perturbations
    
    logpcase      = open(logpcasename,"w+")
    logpcase.write('perturbed buses:\n')
    
    firstline     = lines[0].split()
    
    casename      = firstline[3] + '_line5'
    #casename      = firstline[3] + '_line'
    #casename      = firstline[3] + '_b'
    #casename      = firstline[3] + '_n_0_5'
    #casename      = firstline[3] + '_n_5_5'
    
    firstline[3]  = casename

    newfirstline  = ' '.join(firstline)
    
    newcase.write(newfirstline + '\n')

    while linenum <= len(lines):
        line = lines[linenum-1]
        thisline = line.split()
        if len(thisline) > 0:
          theword = thisline[0]
          if theword[0:4] == 'mpc.':
            newcase.write(line) #writeline
            log.joint("found " + theword + " on line " + str(linenum) + "\n")
            if theword == "mpc.baseMVA":
              boo = thisline[2]
              if boo[len(boo)-1] == ";":
                boo = boo[:len(boo)-1]
              all_data['baseMVA'] = float(boo)
              baseMVA = all_data['baseMVA']
              log.joint(" baseMVA: " + str(baseMVA) + "\n")
            elif theword == 'mpc.bus':
              linenum += 1
              lookingforendofbus = 1
              buses = {}
              numbuses = 0
              slackbus = -1
              IDtoCountmap = {}

              numPload = 0
              sumload = 0
              sumPd = sumQd = 0
              while lookingforendofbus and linenum <= numlines:
                line = lines[linenum-1]
                thisline = line.split()
                if thisline[0] == "];":
                  newcase.write(line) #writeline
                  log.joint("found end of bus section on line " + str(linenum) + "\n")
                  lookingforendofbus = 0
                  break
                numbuses += 1
                if thisline[1] == "3":
                  log.joint(" slack bus: " + thisline[0] + "\n")
                  slackbus = int(thisline[0])
                if thisline[0] != '%':
                  nodeID, nodetype = int(thisline[0]), int(thisline[1])
                  if nodetype != 1 and nodetype != 2 and nodetype != 3 and nodetype != 4:
                    log.joint("bad bus " + thisline[0] + " has type " + thisline[1] + "\n")
                    sys.exit("bad")
                  if nodetype == 4:
                    #log.joint("bus " + thisline[0] + " is isolated\n")
                    numisolated += 1
                    
                  foo = thisline[12]
                  if foo[len(foo)-1] == ';':
                    foo = foo[:len(foo)-1]
                  
                  Pd        = float(thisline[2])
                  PERTURBED = False

                  if ber:
                    if np.random.binomial(size=1,n=1,p=0.3)[0] == 1:
                      PERTURBED = True
                      if randint(0, 1) == 0:
                        newPd = 0.9 * Pd
                      else:
                        newPd = 1.1 * Pd
                      
                  elif normal:
                    PERTURBED = True

                    if normal0:
                      p_mu        = 0
                    elif normal5:
                      p_mu        = 0.05 * Pd
                    
                    p_sigma       = abs(0.05 * Pd)
                    perturbation  = np.random.normal(p_mu,p_sigma,1)[0]
                    newPd         = np.round(Pd + perturbation,2)
                  
                  if PERTURBED:
                    if newPd > 0:
                      log.joint(' bus ' + str(numbuses) + ' nodeID ' + str(nodeID) + ' Pd ' + str(Pd) + ' newPd ' + str(newPd) + '\n')
                      logpcase.write(' bus ' + str(numbuses) + ' nodeID ' + str(nodeID) + ' Pd ' + str(Pd) + ' newPd ' + str(newPd) + '\n')
                      Pd = newPd

                  Vmin  = float(foo)
                  Qd    = float(thisline[3])
                  Gs    = float(thisline[4])
                  Bs    = float(thisline[5])
                  Vbase = float(thisline[9])
                  Vmax  = float(thisline[11])
        
                  thisline[2] = str(Pd)
                  newline = '\t'.join(thisline)
                  newcase.write('\t' + newline + "\n") #writeline

                  buses[numbuses] = Bus(numbuses, nodeID, nodetype, Pd/baseMVA, Qd/baseMVA, Gs/baseMVA, Bs/baseMVA, Vbase,Vmax, Vmin, linenum-1)

                  if nodetype == 1 or nodetype == 2 or nodetype == 3:
                    sumPd += Pd
                    sumQd += Qd

                  IDtoCountmap[nodeID] = numbuses
                  numPload += (Pd > 0)
                else:
                  newcase.write('\t' + line)
                  #newcase.write(line + "\n")
                linenum += 1


              all_data['buses']        = buses
              all_data['numbuses']     = len(buses)
              all_data['sumPd']        = sumPd
              all_data['sumQd']        = sumQd              
              all_data['IDtoCountmap'] = IDtoCountmap
              all_data['slackbus']     = slackbus

              log.joint(" sumloadPd " + str(sumPd) + " numPload " + str(numPload) + "\n")
              log.joint(" sumloadQd " + str(sumQd) + "\n")

              if lookingforendofbus:
                log.stateandquit(" did not find bus data section")
              if slackbus < 0:
                log.joint(" did not find slack bus\n")

              log.joint(" " + str(numbuses)+ " buses\n")
              if numisolated > 0:
                log.joint(" isolated: " + str(numisolated) + "\n")
              
            elif theword == 'mpc.gen':
              linenum += 1
              lookingforendofgen = 1
              gencount = 0
              gens = {}
    
              summaxgenP = summaxgenQ = 0
              while lookingforendofgen and linenum <= numlines:
                line = lines[linenum-1]
                thisline = line.split()
                newcase.write(line) #writeline

                if thisline[0] == "];":
                  all_data['endofgen'] = linenum
                  log.joint(" found end of gen section on line " + str(linenum) + "\n")
                  lookingforendofgen = 0
                  break
                gencount += 1

                #API cases
                if ';' in thisline[9]:
                  thisline[9] = thisline[9].split(';')[0]

                nodeID = int(thisline[0])
                Pg     = float(thisline[1])
                Qg     = float(thisline[2])
                status = int(thisline[7])
                Pmax   = float(thisline[8])
                Pmin   = float(thisline[9])
                Qmax   = float(thisline[3])
                Qmin   = float(thisline[4])

                if status <= 0:
                  status = 0
                else:
                  status = 1

                #log.joint("generator in bus ID " + str(nodeID) )
                if nodeID in IDtoCountmap.keys():
                  idgen = IDtoCountmap[nodeID]
                  gens[gencount] = gen(gencount, nodeID, Pg, Qg, status, Pmax/baseMVA, Pmin/baseMVA, Qmax/baseMVA, Qmin/baseMVA, linenum-1)
                  buses[idgen].addgenerator(log, gencount, gens[gencount])

                  if buses[idgen].nodetype == 2 or buses[idgen].nodetype == 3:  #but not 4
                    summaxgenP += Pmax
                    summaxgenQ += Qmax

                else:
                  log.joint(" generator # " + srt(gencount) + " in nonexistent bus ID " + str(nodeID) + "\n")
                  return 1

                linenum += 1

              if lookingforendofgen:
                log.joint("did not find end of generator section\n")
                return 1
              all_data['gens'] = gens
              all_data['numberofgens'] = len(gens)
              busgencount = 0
              for bus in buses.values():
                busgencount += len(bus.genidsbycount) > 0
              all_data['busgencount'] = busgencount
              log.joint(" number of generators: " + str(all_data['numberofgens']))
              log.joint(" number of buses with gens: " + str(all_data['busgencount']) + "\n")
              all_data['summaxgenP'] = summaxgenP
              all_data['summaxgenQ'] = summaxgenQ              
              log.joint(" summaxPg: " + str(summaxgenP) + " summaxQg: " + str(summaxgenQ) + "\n")

            elif theword == 'mpc.branch':
              lookingforendofbranch = 1
              branchcount = 0
              activebranches = 0
              branches = {}
              linenum += 1
              zerolimit = 0
              while lookingforendofbranch and linenum <= numlines:
                line = lines[linenum-1]
                thisline = line.split()
                
                if thisline[0] == "];":
                  log.joint(" found end of branch section on line " +str(linenum) + "\n")
                  lookingforendofbranch = 0
                  newcase.write(line)
                  break

                branchcount += 1
                f        = int(thisline[0])
                t        = int(thisline[1])
                r        = float(thisline[2])
                x        = float(thisline[3])
                bc       = float(thisline[4])
                rateA    = float(thisline[5])
                rateB    = float(thisline[6])
                rateC    = float(thisline[7])
                ratio    = float(thisline[8])
                angle    = float(thisline[9])
                status   = int(thisline[10])
                minangle = float(thisline[11])

                foo = thisline[12]
                if foo[len(foo)-1] == ';':
                  foo = foo[:len(foo)-1]
        
                maxangle = float(foo)

                if maxangle < minangle:
                  log.stateandquit(" branch # " + str(branchcount) + " has illegal angle constraints\n")

                id_f = IDtoCountmap[f]
                id_t = IDtoCountmap[t]
                
                if pertline:
                  if topbranch == branchcount:
                    status       = 0
                    pertline     = 0
                    list_of_gens_f = buses[id_f].genidsbycount
                    list_of_gens_t = buses[id_t].genidsbycount
                    log.joint(' OFF: branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' genids_f ' + str(list_of_gens_f) + ' genids_t ' + str(list_of_gens_t) + '\n')
                    logpcase.write(' OFF: branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' genids_f ' + str(list_of_gens_f) + ' genids_t ' + str(list_of_gens_t) + '\n')
                    if (buses[id_f].Bs != 0 and topbranch_fdeg == 1) or (buses[id_t].Bs != 0 and topbranch_tdeg == 1):
                      log.joint(' look closely at these cases, \n')
                  
                    thisline[10] = str(status)
                    newline = '\t'.join(thisline)
                    newcase.write('\t' + newline + "\n") #writeline
                  else:
                    newcase.write(line) #writeline
                # if pertline:
                #   if len(buses[id_f].genidsbycount):
                #     status       = 0
                #     pertline     = 0
                #     list_of_gens = buses[id_f].genidsbycount
                #     log.joint(' OFF branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' genids ' + str(list_of_gens) + '\n')
                #     logpcase.write(' OFF branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' genids ' + str(list_of_gens) + '\n')
                #     if buses[id_f].Bs != 0 and len(list_of_gens) == 1:
                #       log.joint(' take care of these cases\n')
                #       breakexit('check!')
                #     breakexit('OFF branch')

                #     thisline[10] = str(status)
                #     newline = '\t'.join(thisline)
                #     newcase.write('\t' + newline + "\n") #writeline

                #   elif len(buses[id_t].genidsbycount):
                #     status       = 0
                #     pertline     = 0
                #     list_of_gens = buses[id_t].genidsbycount
                #     log.joint(' OFF branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' genids ' + str(list_of_gens) + '\n')
                #     logpcase.write(' OFF branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' genids ' + str(list_of_gens) + '\n')
                #     if buses[id_t].Bs != 0 and len(list_of_gens) == 1:
                #       log.joint(' take care of these cases\n')
                #       breakexit('check!')
                #     breakexit('OFF branch')                    
                  
                #     thisline[10] = str(status)
                #     newline = '\t'.join(thisline)
                #     newcase.write('\t' + newline + "\n") #writeline
                #   else:
                #     newline = '\t'.join(thisline)
                #     newcase.write('\t' + newline + "\n") #writeline
                else:
                  newcase.write(line) #writeline

                if status:
                  branches[branchcount] = branch(log,branchcount, f, id_f, t, id_t, r, x, bc, rateA/baseMVA, rateB/baseMVA, rateC/baseMVA, ratio, angle, maxangle, minangle, status, 2*sumPd/baseMVA, linenum-1)
                  zerolimit += (branches[branchcount].constrainedflow == 0)
                  activebranches += 1
        

                  buses[id_f].addfrombranch(log, branchcount)

                  buses[id_t].addtobranch(log, branchcount)

                linenum += 1
              all_data['branches']    = branches
              all_data['numbranches'] = all_data['branchcount'] = branchcount
              log.joint(" branchcount: " + str(branchcount) + " active " + str(activebranches) + "\n")
              log.joint("  " + str(zerolimit) + " unconstrained\n")
            elif theword == 'mpc.gencost':
              lookingforendofgencost = 1
              gencostcount = 1
              linenum += 1
              while lookingforendofgencost and linenum <= numlines:
                line = lines[linenum-1]
                thisline = line.split()
                newcase.write(line) #writeline
                
                if thisline[0] == "];":
                  log.joint(" found end of gencost section on line " +  str(linenum) + "\n")
                  all_data['endofgencost'] = linenum
                  lookingforendofgencost = 0
                  break
                if gencostcount <= gencount:
                  costtype = int(thisline[0])
                  if costtype != 2:
                    log.stateandquit(" cost of generator " + str(gencostcount) + " is not polynomial\n")
                  degree = int(thisline[3]) - 1
                  if degree > 2 or degree < 0:
                    log.stateandquit(" degree of cost function for generator " + str(gencostcount) + " is illegal\n")
                  costvector = [0 for j in range(degree+1)];


                  for j in range(degree+1):
                    boo = thisline[4+j]
                    if boo[len(boo)-1] == ";":
                      boo = boo[:len(boo)-1]
                    costvector[j] = float(boo)

                    costvector[j] *= (baseMVA)**(degree - j)
                      #              print "gen", gencostcount, j, costvector[j]
                  #print costvector

                  gens[ gencostcount ].addcost(log,costvector,linenum)

                else:
                  log.stateandquit(" read " + str(gencostcount) +" gen costs but only " + str(gencount) + " generators\n")
                gencostcount += 1
                linenum += 1
              linenum += 1
          else:
            if 'CASE' in theword:
              newcase.write("%PERTURBED " + casename + "\n")
            else:
              newcase.write(line)
        else:
          newcase.write('\n')
        linenum += 1

    newcase.close()
    logpcase.close()

    return 0

def readvoltsfile(log, all_data):
    voltsfilename = all_data['voltsfilename']
    IDtoCountmap = all_data['IDtoCountmap']
    buses = all_data['buses']
    try:
        f = open(voltsfilename, "r")
        log.joint("reading volts file " + voltsfilename + "\n")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file " + voltsfilename)


    inputvolts = {}
    numread = 0
    
    for linenum in range(len(lines)):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

          if thisline[0] == 'bus':
            angle_rad = float(thisline[5])*math.pi/180
            busid = int(thisline[1])
            inputvolts[busid] = (float(thisline[3]),angle_rad)

            numread += 1
          elif thisline[0] == 'END':
            break
          else:
            print ("illegal input "+ thisline[0] + "\n")
            sys.exit("bye")

    log.joint("read " + str(numread) + " input voltages \n")
    all_data['inputvolts'] = inputvolts
    

def readflowsfile(log, all_data):
    flowsfilename = all_data['flowsfilename']
    IDtoCountmap = all_data['IDtoCountmap']
    buses = all_data['buses']
    try:
        f = open(flowsfilename, "r")
        log.joint("reading flows file " + flowsfilename + "\n")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file " + flowsfilename)

    inputPf = {}
    inputQf = {}    
    inputPt = {}
    inputQt = {}    
    numread = 0

    baseMVA = all_data['baseMVA']
    
    for linenum in range(len(lines)):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

          if thisline[0] == 'branch':
            branchid = int(thisline[1])
            inputPf[branchid] = float(thisline[7])/baseMVA
            inputPt[branchid] = float(thisline[9])/baseMVA
            inputQf[branchid] = float(thisline[11])/baseMVA
            inputQt[branchid] = float(thisline[13])/baseMVA

            numread += 1
          elif thisline[0] == 'END':
            break
          else:
            print ("illegal input " + thisline[0] + " on line " + str(thisline) + "\n")
            sys.exit("bye")

    log.joint("read " + str(numread) + " input flows\n")

    all_data['inputPf'] = inputPf
    all_data['inputPt'] = inputPt
    all_data['inputQf'] = inputQf
    all_data['inputQt'] = inputQt
    
    
def writegv(log, all_data, gvfilename):
  try:
    f = open(gvfilename, "w")
    log.joint("writing to gv file " + gvfilename + "\n")
  except:
    log.stateandquit("cannot open file " + gvfilename)
  f.write("graph {\n")

  for bus in all_data['buses'].values():
    f.write("     " + str(bus.nodeID)+";\n")
  for branch in all_data['branches'].values():
    f.write("     " + str(branch.f)+" -- " + str(branch.t)+";\n")

  f.write("}\n")
  f.close()
  sys.exit()


def generateinputcs(log,all_data):
  log.joint("  generating input c,s values\n")
  inputvolts = all_data['inputvolts']
  inputcc = {}
  inputcs = {}

  buses = all_data['buses']
  branches = all_data['branches']
  IDtoCountmap = all_data['IDtoCountmap']
  for busid in inputvolts:
    M = inputvolts[busid][0]
    #    log.joint(str(busid) +" is in input volts with M " +  str(M) + "\n")
    inputcc[busid] = M*M

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    if f in inputvolts and t in inputvolts:
      Mf = inputvolts[f][0]
      Mt = inputvolts[t][0]
      af = inputvolts[f][1]
      at = inputvolts[t][1]
      angle = af - at
      inputcs[f,t] = (Mf*Mt*math.cos(angle),Mf*Mt*math.sin(angle))
      inputcs[t,f] = (Mf*Mt*math.cos(angle),-Mf*Mt*math.sin(angle))


  all_data['inputcc'] = inputcc 
  all_data['inputcs'] = inputcs 
  


def generateinputeandf(log,all_data):
  log.joint("  generating input e,f values\n")
  inputvolts = all_data['inputvolts']
  inputve = {}
  inputvf = {}

  buses = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap']
  for busid in inputvolts:
    M = inputvolts[busid][0]
    A = inputvolts[busid][1]
    #    log.joint(str(busid) +" is in input volts with M " +  str(M) + "\n")
    inputve[busid] = M*math.cos(A)
    inputvf[busid] = M*math.sin(A)


  all_data['inputve'] = inputve
  all_data['inputvf'] = inputvf
  
def readdigits(log, all_data):
  Lfilename = all_data['Lfilename']
  log.joint("reading L file " + Lfilename + "\n")

  buses = all_data['buses']

  try:
    f = open(Lfilename, "r")
    lines = f.readlines()
    f.close()
  except:
    log.stateandquit("cannot open file " + Lfilename)
    sys.exit("failure")

  L = {}
  IDtoCountmap = all_data['IDtoCountmap'] 

  for bus in buses.values():
    L[bus] = 0
  for linenum in range(len(lines)):
    line = lines[linenum].split()
    if line[0] == 'default':
      lvalue = int(line[1])
      log.joint(" default L: " + str(lvalue) + "\n")
      for bus in buses.values():
        L[bus] = lvalue
    elif line[0] == 'END':
      break
    else:
      ind = int(line[0])
      L[buses[ IDtoCountmap[ind] ] ]= int(line[1])

  all_data['L'] = L

    
