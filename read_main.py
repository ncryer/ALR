import numpy as np

"""
Read ALR files and update a contig data structure
"""


    
def scan_file(filename, ContigHolder ):
    
    """
    Use buffered I/O to loop over the ALR file and parse each line
    
    Update an aggregator object as we go along
    """
    
    with open(filename) as f:
        
        for line in f:
            
            ContigHolder.build_contig( parse_line( line ) )
    
    ContigHolder.current_contig.summarize_counts()
    ContigHolder.contigs.append(ContigHolder.current_contig)
    ContigHolder.current_contig = None
            
    return ContigHolder
            
            
def parse_line( line ):
    
    line = line.strip().split( "\t" )
    
    # Determine line type

    
    if line[0].startswith(">"):
        
        # New contig! End the old one 
        
        return ("name", line[0])
    
    if line[0] == "maj":
        
        # Header information
        
        return ("header", line)
    
    if line[1] == "M":
        
        # Standard locus
        
        return ("maj", line)
    
    if line[1] == "P":
        
        # Polymorphic site
        
        return ("poly", line)
        
    
    
class ContigStruct(object):
    
    """
    Aggregates contig objects
    """    
    
    def __init__(self):
        
        self.organism_name = ""
        
        self.contigs = []
        
        self.current_contig = False
        
        self.header = False
        
        self.n_columns = 0
        
    def parse_header(self, line):
        
        inf = line[2:]
        
        self.n_columns = len(inf)
        
        self.organism_name = inf[0].split("|")[0]
        
        self.header = [item.split("|")[1] for item in inf]
    
    def build_contig(self, linetup):
        """
        Accepts a tuple where indx 0 specifies the type of line
        """
        
        linetype = linetup[0]
        line = linetup[1]
        
        if linetype == "name": 
            # New contig starts, finish old one if necessary
            
            if self.current_contig:
                
                # Finish this and add it to the contig list
                
                self.current_contig.summarize_counts()
                
                self.contigs.append( self.current_contig )
                
                self.current_contig = None
                
                self.current_contig = Contig( int(line[7:]), self.n_columns )
                
            else:
                
                # First contig
                # TODO: Clean this up
                self.current_contig = Contig( int(line[7:]), 18 )
        
        if linetype == "header" and self.header == False:
            
            # Add header information, otherwise let it pass
            
            self.parse_header(line)
        
        if linetype == "maj":
            
            self.current_contig.add_majline(line)
        
        if linetype == "poly": 
            
            self.current_contig.add_polyline(line)
    
    def write_count_matrix(self):
        
        import csv 
        
        f = open("countmat.csv", "w")
        
        outheader = ["Name"] + self.header
        
        writer = csv.writer(f, dialect="excel")
        
        writer.writerow(outheader)

        for contig in self.contigs:
            
            writer.writerow( [contig] + list(contig.summarized_counts) )
        
        f.close()
        
        print "Done."
            
            

class Contig:
    
    def __init__(self, name, n_columns):
        
        self.name = name
        self.counts = np.asarray([0]*n_columns)
        self.seq = ""
        self.seqlen = 0
        self.summarized_counts = self.counts
        self.polysites = []
        self.nucdict = {0 : "A", 1: "C", 2 : "G", 3: "T"}
    
    def add_majline(self, line):
        
        """ Extract information from an M locus"""
        
        # Take care of variable length sequence
        
        if len(line[2:]) < len(self.counts):
            
            dif = len(self.counts) - len(line[2:])
            
            line = line + [0]*dif
        
        # Add nucleotide to sequence
        nuc = str(line[0])
        self.seq += nuc
        # Increment seq length
        self.seqlen += 1
        
        # Convert the reads to integer type
        reads = np.array( map(int, line[2:]) )

        
        # Increment the counts

        self.counts += reads

                

    def add_polyline(self, line):

        """ Extract information from a polymorphic site """
        
        # Take care of variable length sequence
        
        if len(line[2:]) < len(self.counts):
            
            dif = len(self.counts) - len(line[2:])
            
            for i in range(dif):
                
                line.append("0[0/0/0/0]")
        
        # Add nucleotide to sequence
        nuc = str(line[0])
        self.seq += nuc
        self.seqlen += 1
        
        # TODO: Find the minority allele 
        
        # Get majority reads
        inp = line[2:]
        
        reads = np.array(map(int, [ inp[i][0] for i in range(len(inp)) ] )) 
        self.counts += reads
        
    def summarize_counts(self):
        
        """
        Average counts by sequence length in each column
        """
        
        self.summarized_counts = np.round( self.counts / np.float(self.seqlen) )
        

    def print_info(self):
        # For debugging   
        print "NAME: " + str(self.name) + "\n"
        print "COUNTS: " + str(self.counts) + "\n"
        print "SUMMARIZED COUNTS: " + str(self.summarized_counts) + "\n"
        print "SEQ: " + str(self.seq) + "\n"
        print "SEQLEN: " + str(self.seqlen) + "\n"
        
    def __repr__(self):
        
        return "Contig" + str(self.name)        
        
if __name__ == "__main__":

    ContigHolder = ContigStruct()
    
    scan_file('Galaxy10-[F__UTR5_region_(ALR_format)].alr', ContigHolder )
    
    #ContigHolder.write_count_matrix()
