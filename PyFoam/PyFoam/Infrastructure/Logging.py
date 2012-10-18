#  ICE Revision: $Id: /local/openfoam/Python/PyFoam/PyFoam/Infrastructure/Logging.py 1906 2007-08-28T16:16:19.392553Z bgschaid  $ 
"""Writes Logfiles"""

try:
    import logging
    hasLogging=True
except ImportError:
    # For Python-versions < 2.3
    print "Warning: old python-version. No logging-support"
    hasLogging=False
    
from Hardcoded import assertDirectory,logDirectory
from os import path,uname

from PyFoam import configuration as config

_definedLoggers=[]

def _getLoggingLevel(name):
    """Gets the logging level value from its name"""
    level=config().get("Logging","default")

    try:
        level=config().get("Logging",name)
    except:
        pass

    value=logging.INFO
    try:
        value=getattr(logging,level)
    except AttributeError,reason:
        print "WARNING: Wrong specification of debug level "+level+" for log "+name
        
    return value

class DummyLogger:
    def __init__(self):
        pass

    def info(self,arg):
        pass
    
def foamLogger(name="general"):
    """
    @param name: name of the logfile
    @return: a logger that is correctly set up for pyFoam
    """

    if not hasLogging:
        return DummyLogger()
    
    log=logging.getLogger(name)
    if not (name in _definedLoggers):
        assertDirectory(logDirectory())
        lname=path.join(logDirectory(),name)
        #        rot=logging.TimedRotatingFileHandler(lname,when="M",interval=2,backupCount=5)
        rot=logging.FileHandler(lname)
        machine=uname()[1].split(".")[0]
        rot.setFormatter(logging.Formatter(fmt="%(asctime)s "+("%15s" % machine)+":%(process)-6d %(levelname)-8s %(message)s - in %(filename)s:%(lineno)d"))
        log.addHandler(rot)
        log.setLevel(_getLoggingLevel(name))
        _definedLoggers.append(name)
        
    return log
