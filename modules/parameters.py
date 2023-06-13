


from dataclasses import dataclass


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



