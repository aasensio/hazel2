class Error(Exception):
   """Base class for other exceptions"""
   pass

class NumericalErrorSIR(Error):
   """Raised when there was a problem with SIR"""
   pass

class NumericalErrorHazel(Error):
   """Raised when there was a problem with Hazel"""
   pass