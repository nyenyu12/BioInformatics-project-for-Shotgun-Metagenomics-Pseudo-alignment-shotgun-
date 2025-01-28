import re

class NoRecordsInDataFile(Exception):    
    def __init__(self, message, errors):            
        super().__init__(message)

class InvalidRecordData(Exception):    
    def __init__(self, message, errors):            
        super().__init__(message)

class Record(object):
       
    def __init__(self, sections):
        if (len(sections) == 0) or (len(sections[0]) == 1):
            raise InvalidRecordData("The data given to construct record has invalid dimensions.")
        
        self.identifier = sections[0][0]
        self.__sections = {}
        
        for section_header, section_data in sections:
            if section_header in self.__sections:
                raise InvalidRecordData(f"Section header: {section_header} has appeared twice in the given data.")
            
            self.__sections[section_header] = section_data
    
        
class RecordContainer(object):
    RECORD_SECTIONS = {}


    def get_record_re_pattern():
        re_pattern = []
        
        for section_header, section_legal_chars in Record.RECORD_SECTIONS.keys(), Record.RECORD_SECTIONS.values():
            re_pattern.append(f'^{section_header}')
            if len(section_header) == 0:
                re_pattern.append(f'+{section_legal_chars}')
            else:
                re_pattern.append(f'*{section_legal_chars}')
            
        return re_pattern


class DataFile(object):
    
    EXTENSIONS = ()

    def __init__(self, data_file_path):
        self.__records = {}
        self.parse_file(data_file_path)

    def parse_file(self, data_file_path):
        with open(data_file_path, 'rb') as data_file:
            self.parse_records(data_file)
            
        if len(self.__records) == 0:
            raise NoRecordsInDataFile
        
    @property    
    def extensions():
        return DataFile.EXTENSIONS
    
    def find_header(data_file):
        raise NotImplementedError
    
    def find_record_section_headers(data_file):
        pass
    
    def parse_records(self, data_file):
        re.finditer()
        
        raise NotImplementedError
    
        
class FastaFile(DataFile):
    

    EXTENSIONS = set('.fa', '.fa.gz')
    CHARS_TO_IGNORE = set(whitespace)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)