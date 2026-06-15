class FudanCourse():
    def __init__(self,course_id,place) -> None:
        self.id=course_id
        self.place=place
        self.PNP=True # 自动初始化的属性
    def info(self):
        print("ID: "+self.id+"; Place: "+self.place+"; P/NP:"+str(self.PNP))
    def change_place(self,place):
        self.place=place
        
class Teacher:
    def __init__(self,department,name) -> None:
        self.department=department  
        self.name=name
    def change_name(self,name):
        self.name=name
    def change_department(self,dpt):
        self.department=dpt

class PhysicsCourse(FudanCourse):
    def __init__(self, course_id, place) -> None:
        super().__init__(course_id, place) # 让实例包含父类的所有属性
        self.experiment=True # 新的属性
        self.teacher=Teacher('Physics','?') # 用类的实例作为属性
    def info(self):
        print("ID: "+self.id+"; Place: "+self.place+"; P/NP:"+str(self.PNP)+
              "; Experiments:"+str(self.experiment)+"; Teacher: "+self.teacher.name+
              " from "+self.teacher.department)