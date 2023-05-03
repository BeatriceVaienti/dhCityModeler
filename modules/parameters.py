


from dataclasses import dataclass


"""
The best case is when we know the number of floors and the floor height. 
in this case the floor height is not used and is merely inferred for later used

if we know the height but not the number of floors or the floor height, 
we can infer the number of floors by a default floor height and readjust the floor height to obtain an integer number of floors

if we know the number of floors but not the height or the floor height, 
we can infer the height by multiplying the number of floors by a default floor height

if we know the number of floors and the floor height then the overall height is simply inferrecd 
by multiplying the number of floors by the floor height

if we don't know anything, the number of floor is assigned randomly in a given range and the height is inferred by multiplying the number of floors by a default floor height 
"""

@dataclass
class Source:  
    name: str = None
    type: str = None
    notes: str = None

    def create_dictionary(self):
        return {"name":self.name, "type":self.type, "notes":self.notes}

    def fill(self, field, value):
        if field == "name":
            self.name = value
        if field == "type":
            self.type = value
        if field == "notes":
            self.notes = value

@dataclass
class Paradata:
    author: str = None
    date: str = None
    comments: str  = None
    uncertainty: str = None
    version: str = 0

    def create_dictionary(self):
        return  {"author":self.author, "date":self.date, "comments":self.comments, "uncertainty":self.uncertainty, "version":self.version}

    def fill(self, field, value):
        if field == "author":
            self.author = value
        if field == "date":
            self.date = value
        if field == "comments":
            self.comments = value
        if field == "uncertainty":
            self.uncertainty = value 
        if field == "version":
            self.version = value




@dataclass
class Parameter:  
    name: str = None
    value: str = None
    sources: list = None
    paradata: str = None
    
    def create_dictionary(self):
        return {"value":self.value, "sources":self.sources, "paradata": self.paradata } 
    
    def fill(self, field, value):
        if field == "value":
            self.value = value
        if field == "sources":
            self.sources = value
        if field == "paradata":
            self.paradata = value



