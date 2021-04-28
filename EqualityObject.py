class EqualityObject:
    def __init__(self,node,length,equalitycounter,supportlist):
        self.node=node
        self.len=length
        self.eq=equalitycounter
        self.support=[]
        self.support.append(supportlist)

    def add_supported_read(self,supported_read):
        self.support.append(supported_read)
    def get_node(self):
        return self.node
    def get_length(self):
        return self.len
    def get_support(self):
        return self.support
    def get_equality(self):
        return self.eq
    def increment_eq(self):
        self.eq += 1
    def __str__(self):
        s="".join((str(e)+" ") for e in self.support)
        return str(self.len)+", "+str(self.eq)+", "+s+"\n"