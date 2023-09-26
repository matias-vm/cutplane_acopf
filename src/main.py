#!/usr/bin/python

import sys
import os
import numpy as np
import time
import reader
from myutils import breakexit
from versioner import *
from log import danoLogger
from cutplane import gocutplane

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    casefilename                 = 'NONE'
    lpfilename                   = ""
    lpfilename_cuts              = ""
    linenum                      = 0

    
    jabrcuts                     = 0
    i2cuts                       = 0
    losscuts                     = 0
    limitcuts                    = 0
    objective_cuts               = 0    
    loud_cuts                    = 0


    jabr_inequalities            = 0
    i2_inequalities              = 0
    limit_inequalities           = 0
    loss_inequalities           = 0



    threshold                    = 1e-5
    threshold_i2                 = 1e-1
    tolerance                    = 1e-5
    threshold_objcuts            = 1e-5
    threshold_dotprod            = 5e-1
    most_violated_fraction       = 1
    dropjabrs                    = 0
    droploss                     = 0
    dropi2                       = 0
    cut_age_limit                = 20
    solver_method                = 2
    primal_bound                 = 'NONE'
    crossover                    = 0
    timelimit                    = 'NONE'
    max_rounds                   = 100
    cut_analysis                 = 0

    most_violated_fraction_i2    = 1
    i2_def_threshold             = 1e-3
    mincut                       = 0
    mincut_active                = 0
    mincut_reactive              = 0
    mincut_switch                = 0

    dographics                   = 0
    
    obbt                         = 0

    hybrid                       = 0

    jabr_validity                = 0
    i2_validity                  = 0
    loss_validity                = 0
    limit_validity               = 0

    getsol                       = 0

    fixflows                     = 0
    fixcs                        = 0
    tol_fix                      = 1e-5
    knitro_sol                   = 0
    matpower_sol                 = 0
    writeACsol                   = 0

    FeasibilityTol               = 1e-5


    linear_objective             = 0

    writecuts                    = 0
    addcuts                      = 0


    droplimit                    = 0
    most_violated_fraction_limit = 1
    threshold_limit              = 1e-5


    writelps                     = 0

    ftol                         = 1e-3
    ftol_iterates                = 5
    max_time                     = 200


    

    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename = thisline[1]

            elif thisline[0] == 'lpfilename':
                lpfilename = thisline[1]

            elif thisline[0] == 'lpfilename_cuts':
                lpfilename_cuts = thisline[1]
            
            elif thisline[0] == 'jabr_inequalities':
                jabr_inequalities = 1

            elif thisline[0] == 'cut_analysis':
                cut_analysis = 1
                
            elif thisline[0] == 'dropjabrs':
                dropjabrs = 1

            elif thisline[0] == 'droploss':
                droploss = 1

            elif thisline[0] == 'dropi2':
                dropi2 = 1                

            elif thisline[0] == 'limit_inequalities':
                limit_inequalities = 1

            elif thisline[0] == 'solver_method':
                solver_method = int(thisline[1])

            elif thisline[0] == 'objective_cuts':
                objective_cuts = 1

            elif thisline[0] == 'linear_objective':
                linear_objective = 1

            elif thisline[0] == 'i2cuts':
                i2cuts           = 1

            elif thisline[0] == 'i2_def_threshold':
                i2_def_threshold = float(thisline[1])

            elif thisline[0] == 'most_violated_fraction_i2':
                most_violated_fraction_i2 = float(thisline[1])                
                
            elif thisline[0] == 'cut_age_limit':
                cut_age_limit = int(thisline[1])

            elif thisline[0] == 'threshold':
                threshold = float(thisline[1])

            elif thisline[0] == 'threshold_i2':
                threshold_i2 = float(thisline[1])
     
            elif thisline[0] == 'threshold_objcuts':
                threshold_objcuts = float(thisline[1])

            elif thisline[0] == 'threshold_dotprod':
                threshold_dotprod = float(thisline[1])
                
            elif thisline[0] == 'tolerance':
                tolerance = float(thisline[1])

            elif thisline[0] == 'mincut':
                mincut          = 1
                mincut_active   = 1

            elif thisline[0] == 'mincut_reactive':
                mincut          = 1
                mincut_active   = 0
                mincut_reactive = 1

            elif thisline[0] == 'mincut_switch':
                mincut_switch   = 1

            elif thisline[0] == 'most_violated_fraction':
                most_violated_fraction = float(thisline[1])

            elif thisline[0] == 'primal_bound':
                primal_bound = float(thisline[1])

            elif thisline[0] == 'crossover':
                crossover = int(thisline[1])

            elif thisline[0] == 'max_time':
                max_time  = float(thisline[1])

            elif thisline[0] == 'max_rounds':
                max_rounds = int(thisline[1])

            elif thisline[0] == 'dographics':
                dographics = 1

            elif thisline[0] == 'obbt':
                obbt = 1

            elif thisline[0] == 'hybrid':
                hybrid = 1

            elif thisline[0] == 'i2_inequalities':
                i2_inequalities = 1

            elif thisline[0] == 'jabr_validity':
                jabr_validity = 1

            elif thisline[0] == 'i2_validity':
                i2_validity = 1

            elif thisline[0] == 'loss_validity':
                loss_validity = 1

            elif thisline[0] == 'fixflows':
                fixflows = 1

            elif thisline[0] == 'knitro_sol':
                knitro_sol = 1

            elif thisline[0] == 'fixcs':
                fixcs    = 1

            elif thisline[0] == 'tol_fix':
                tol_fix    = float(thisline[1])

            elif thisline[0] == 'getsol':
                getsol     = 1

            elif thisline[0] == 'loss_inequalities':
                loss_inequalities  = 1

            elif thisline[0] == 'most_violated_fraction_loss':
                most_violated_fraction_loss  = float(thisline[1])

            elif thisline[0] == 'FeasibilityTol':
                FeasibilityTol  = float(thisline[1])

            elif thisline[0] == 'writecuts':
                writecuts  = 1

            elif thisline[0] == 'addcuts':
                addcuts    = 1

            elif thisline[0] == 'droplimit':
                droplimit  = 1

            elif thisline[0] == 'most_violated_fraction_limit':
                most_violated_fraction_limit    = 1

            elif thisline[0] == 'threshold_limit':
                threshold_limit     = float(threshold_limit)

            elif thisline[0] == 'limit_validity':
                limit_validity  = 1

            elif thisline[0] == 'limitcuts':
                limitcuts    = 1

            elif thisline[0] == 'writelps':
                writelps     = 1

            elif thisline[0] == 'ftol':
                ftol          = float(thisline[1])

            elif thisline[0] == 'ftol_iterates':
                ftol_iterates = int(thisline[1])

            elif thisline[0] == 'matpower_sol':
                matpower_sol  = 1

            elif thisline[0] == 'knitro_sol':
                knitro_sol    = 1

            elif thisline[0] == 'writeACsol':
                writeACsol    = 1

            elif thisline[0] == 'jabrcuts':
                jabrcuts      = 1

            elif thisline[0] == 'losscuts':
                losscuts      = 1

            elif thisline[0] == 'loud_cuts':
                loud_cuts     = 1

            elif thisline[0] == 'END':
                break
                
            else:
                sys.exit("illegal input " + thisline[0] + "bye")


        linenum += 1

    all_data                    = {}
    all_data['casefilename']    = casefilename
    all_data['casename']        = casefilename.split('data/')[1].split('.m')[0] 
    all_data['lpfilename']      = all_data['casename'] + '.lp'
    all_data['lpfilename_cuts'] = all_data['casename'] + '_cuts.lp'    

    if len(lpfilename):
        all_data['lpfilename']  = lpfilename
    
    if len(lpfilename_cuts):
        all_data['lpfilename_cuts'] = lpfilename_cuts
    

    all_data['jabr_inequalities']           = jabr_inequalities
    all_data['limit_inequalities']          = limit_inequalities
    all_data['i2_inequalities']             = i2_inequalities


    all_data['max_rounds']      = max_rounds
    all_data['solver_method']   = solver_method
    all_data['cut_analysis']    = cut_analysis
    
    #cuts
    all_data['initial_threshold']      = threshold
    all_data['threshold']              = threshold
    all_data['tolerance']              = tolerance
    
    #parallel cuts
    all_data['threshold_dotprod'] = threshold_dotprod
    
    all_data['tight_cuts'] = {}
    all_data['tight_cuts_fraction'] = {}
    
    all_data['cut_age_limit'] = cut_age_limit
    all_data['droploss'] = droploss

    
    all_data['linear_objective']   = linear_objective
    all_data['hybrid']             = hybrid

    if linear_objective or hybrid:
        all_data['objective_cuts']     = objective_cuts
        all_data['obj_cuts']           = {}
        all_data['num_objective_cuts'] = 0
        all_data['threshold_objcuts']  = threshold_objcuts

    all_data['jabrcuts'] = jabrcuts
    if jabrcuts:
        all_data['most_violated_fraction'] = most_violated_fraction
        all_data['ID_jabr_cuts']           = 0
        all_data['jabr_cuts']              = {}
        all_data['num_jabr_cuts']          = 0
        all_data['num_jabr_cuts_rnd']      = {}
        all_data['num_jabr_cuts_added']    = 0
        all_data['num_jabr_cuts_dropped']  = 0
        all_data['dropped_jabrs']          = []    
        all_data['dropjabrs']              = dropjabrs
        all_data['jabr_cuts_info']         = {}
        all_data['jabr_cuts_info_updated'] = {}
        all_data['max_error']              = 0

    all_data['jabr_validity']          = jabr_validity

    all_data['limitcuts'] = limitcuts
    if limitcuts:
        all_data['most_violated_fraction_limit']    = most_violated_fraction_limit
        all_data['ID_limit_cuts']                   = 0
        all_data['limit_cuts']                      = {}
        all_data['num_limit_cuts']                  = 0
        all_data['num_limit_cuts_rnd']              = {}
        all_data['num_limit_cuts_added']            = 0
        all_data['num_limit_cuts_dropped']          = 0
        all_data['droplimit']                       = droplimit
        all_data['limit_cuts_info']                 = {}
        all_data['limit_cuts_info_updated']         = {}
        all_data['threshold_limit']                 = threshold_limit
        all_data['dropped_limit']                   = []
    
        all_data['max_error_limit']                 = 0
    
    all_data['limit_validity']                  = limit_validity    

    all_data['i2cuts'] = i2cuts
    if i2cuts:
        all_data['most_violated_fraction_i2'] = most_violated_fraction_i2
        all_data['ID_i2_cuts']                = 0
        all_data['i2_cuts']                   = {}
        all_data['num_i2_cuts']               = 0
        all_data['num_i2_cuts_rnd']           = {}
        all_data['num_i2_cuts_added']         = 0
        all_data['num_i2_cuts_dropped']       = 0
        all_data['dropi2']                    = dropi2
        all_data['i2_cuts_info']              = {}
        all_data['i2_cuts_info_updated']      = {}
        all_data['threshold_i2']              = threshold_i2
        all_data['dropped_i2']                = []
        all_data['i2_def_threshold']          = i2_def_threshold

        all_data['max_error_i2']              = 0

    all_data['i2_validity']               = i2_validity
        
    all_data['primal_bound']                  = primal_bound
    all_data['crossover']                     = crossover
    all_data['max_time']                      = max_time

    all_data['mincut']                        = mincut
    all_data['mincut_active']                 = mincut_active
    all_data['mincut_reactive']               = mincut_reactive
    all_data['mincut_switch']                 = mincut_switch

    all_data['dographics']                    = dographics
    
    if all_data['dographics']:
        all_data['coordsfilename'] = None

    all_data['obbt']                          = obbt
    if all_data['obbt']:
        all_data['obbt_itcount'] = 0
        all_data['mincut_jabrs'] = []
        all_data['mincut_i2']    = []


    all_data['FeasibilityTol']                = FeasibilityTol
    all_data['fixflows']                      = fixflows
    all_data['fixcs']                         = fixcs
    all_data['tol_fix']                       = tol_fix
    all_data['matpower_sol']                  = matpower_sol
    all_data['knitro_sol']                    = knitro_sol
    all_data['writeACsol']                    = writeACsol



    all_data['loss_inequalities'] = loss_inequalities
    all_data['losscuts']          = losscuts

    if loss_inequalities or losscuts:
        all_data['loss_cuts']                   = {}
        all_data['num_loss_cuts']               = 0
        all_data['num_loss_cuts_added']         = 0
        all_data['num_loss_cuts_dropped']       = 0
        all_data['dropped_loss']                = []
        all_data['most_violated_fraction_loss'] = most_violated_fraction_loss

        all_data['max_error_loss']              = 0

    all_data['loss_validity']               = loss_validity


    all_data['writecuts']                     = writecuts
    all_data['addcuts']                       = addcuts
    all_data['writelps']                      = writelps

    all_data['ftol']                          = ftol
    all_data['ftol_iterates']                 = ftol_iterates

    all_data['loud_cuts']                     = loud_cuts


    return all_data
        

if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) < 2:
        print ('Usage: main.py file.config [logfile]\n')
        exit(0)

    T0 = time.time()
    mylogfile = "main.log"

    if len(sys.argv) == 3:
        #mylogfile = 'logfiles/' + sys.argv[2]
        mylogfile = sys.argv[2] #trial.sh

    log = danoLogger(mylogfile)
    stateversion(log)

    all_data = read_config(log,sys.argv[1])    
    all_data['T0'] = T0

    readcode = reader.readcase(log,all_data,all_data['casefilename'])
    #readcode = newreader.readcase(log,all_data,all_data['casefilename'])    
    #readcode = pglibreader.readcase(log,all_data,all_data['casefilename'])
    #readcode = reader_pertline.readcase(log,all_data,all_data['casefilename'])
    #breakexit('next instance')
    #exit(0)

    gocutplane(log,all_data)
    
    log.closelog()
    
