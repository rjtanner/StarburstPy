# -*- coding: utf-8 -*-
"""
    A method for managing output and error messages from StarburstPy.

    author: Ryan Tanner
    email: ryan.tanner@nasa.gov

"""

"""
    For future development:
        
        Impliment a rank = -1 where errors are suppressed (not raised) and 
        the code continues to run. Useful when running multiple models.
    
        Have the buffer written to a file.

"""


class sb_messages(object):
    """
    Catches and controls all errors, warnings, and messages.
    Outputs are printed based on the rank of the verbosity.
    Regardless of rank all output is saved in a buffer, and can be 
    displayed after code terminates.
    
    rank = 0 : No output
    rank = 1 : Only errors
    rank = 2 : Errors and warnings
    rank = 3 : Errors, warnings, and messages
    """
    
    def __init__(self, rank = None, function = None, write_to = None, stay = True):
        
        self.buffer = {}
        self.rank = rank if rank is not None else 2
        self.erase_buffer()
        
        
    def erase_buffer(self):
        
        self.buffer['errors'] = []
        self.buffer['warnings'] = []
        self.buffer['messages'] = []
        
        
    def message(self, text, function = None):
        """
        Print a general output message and add it to the buffer.
        """
        if self.rank > 2:
            self.print_single(' : ', text, function=function)
        self.buffer['messages'].append('message: {0}: {1}'.format(function, text))


    def warning(self, text, function = None):
        """
        Print a warning and add it to the buffer.
        """
        if self.rank > 1:
            self.print_single('warning', text, function=function)
        self.buffer['warning'].append('warning: {0}: {1}'.format(function, text))
        
        
    def error(self, text, function = None, exception = None):
        """
        Print an error and add it to the buffer. Raise the exception.
        """
        if self.rank > 0:
            self.print_single('error', text, function=function)
        self.buffer['errors'].append('error: {0}: {1}'.format(function, text))
        if exception is None:
                exception = SystemExit
        raise exception(text)
        
        
    def print_single(self, out_type, text, function=None):
        """
        Function to print a single message.
        """
        if function is None:
            function = self.function
        print('{0} in {1}: {2}'.format(out_type, function, text))
        

    def print_buffer(self, erase = True):
        """
        Prints all messages and erases the buffer.
        """
        self.print_messages()
        self.print_warnings()
        self.print_errors()
        if erase:
            self.erase_buffer()
        
        
    def print_errors(self):
        """
        Print out all the errors.
        """
        for text in self.buffer['errors']:
            print(text)
        
    
    def print_warnings(self):
        """
        Print out all the warnings.
        """
        for text in self.buffer['warnings']:
            print(text)
            
            
    def print_messages(self):
        """
        Print out all the messages.
        """
        for text in self.buffer['messages']:
            print(text)
    
    
        
        