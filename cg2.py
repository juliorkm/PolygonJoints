# Second Computer Graphics assignment 2017.2
# @author Julio Rama
# @date 05/11/2017
#
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from random import *
import numpy as np
import geometry
import matrix

# Global variable that stores the current height of the window.
# Required for using the y coordinate of the mouse.
global windowHeight
windowHeight = 600
# Global variable that stores the current width of the window.
global windowWidth
windowWidth = 800

# Boolean that says whether the help menu is active or not.
global helpActive
helpActive = False

# Float that dictates whether the polygons are translucent or not
global alpha
alpha = 1

# Boolean that dictates whether the polygons have a stroke or not
global stroke
stroke = False
# Int that dictates the size of the polygons' stroke
global strokeSize
strokeSize = 2

# Boolean that dictates whether the polygons are in grayscale or not
global grayScale
grayScale = False

# Boolean that dictates whether the polygons' triangulations are explicitly drawn or not
global showTriangles
showTriangles = False

# Boolean that dictates whether the nails will be drawn or not
global showNails
showNails = True

# Boolean and coordinates for drawing temporary lines
isDrawingPolygon = False
global mouseX
global mouseY

# Global list that stores polygons drawn by the user.
global polygons
polygons = []
# Global list that stores nails placed by the user.
global nails
nails = []

# Point in which a polygon was grabbed.
global grabbedPoint
grabbedPoint = None
# Polygon that is being held by the user.
global grabbedPolygon
grabbedPolygon = None
# Point which a polygon and its children should rotate around.
global rotationPoint
rotationPoint = None

# List that stores all intersection points' coordinates.
# Since the intersections cannot change position (lines don't move or cease to exist),
# we can keep storing all new intersections as we insert new lines.
global lineIntersections
lineIntersections = []


# Instead of using the geometry.py's Point class, I'm overriding the comparichild
# method, as to have a margin for error.
class Vertex(geometry.Point):
	def __eq__(self, other):
		return isinstance(other, self.__class__) and \
			   abs(self.x - other.x) < 5 and \
			   abs(self.y - other.y) < 5
			   

# Class that represents the polygons used in the tree. It is also used for representing the polygons
# in a way that OpenGL can draw them.
class Polygon(object):
	# A polygon is initiated by giving it its first vertex (from the first mouse click).
	# It also acquires its own randomly generated color, used when drawing it.
	def __init__(self, x, y):
		self.vertexes = [Vertex(x,y)]
		# Sets a random color to the polygon.
		self.randomizeColor()

		self.father = None # every polygon is initially a root
		self.fatherNail = None # if a polygon has a father, it also stores what nail attached it to him
		self.children = [] # polygons can have multiple children
		self.triangles = [] # polygons have triangles that will be used to draw them
		
	# Function that gives the polygon a random color that isn't too dark or too bright.
	def randomizeColor(self):
		r = random()
		g = random()
		b = random()
		# If the randomly generated color is too dark, it is brightened.
		if r < .5 and g < .5 and b < .5:
			aux = 1.5 - (r+g+b)
			r = r + aux/3
			g = g + aux/3
			b = b + aux/3
		# If the randomly generated color is too bright, it is darkened.
		elif r > .8 and g > .8 and b > .8:
			aux = (r+g+b) - 2.5
			r = r - aux/3
			g = g - aux/3
			b = b - aux/3

		# Cap all values to 1.
		if r > 1: r = 1
		if g > 1: g = 1
		if b > 1: b = 1

		self.color = (r,g,b)
	
	# Triangulation method
	# Uses the GLU tesselator to create a list of triangles based on the polygon's vertexes.
	# These triangles are represented as lists of x,y tuples that are used when rendering the polygon.
	def triangulate(self):
		vertexes = []
		def edgeFlagCallback(param1, param2): pass
		def beginCallback(param=None):
			vertexes = []
		def vertexCallback(vertex, otherData=None):
			vertexes.append(vertex[:2])
		def combineCallback(vertex, neighbors, neighborWeights, out=None):
			out = vertex
			return out
		def endCallback(data=None): pass
		tess = gluNewTess()
		gluTessProperty(tess, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_ODD)
		gluTessCallback(tess, GLU_TESS_EDGE_FLAG_DATA, edgeFlagCallback) # forces triangulation of polygons (i.e. GL_TRIANGLES) rather than returning triangle fans or strips
		gluTessCallback(tess, GLU_TESS_BEGIN, beginCallback)
		gluTessCallback(tess, GLU_TESS_VERTEX, vertexCallback)
		gluTessCallback(tess, GLU_TESS_COMBINE, combineCallback)
		gluTessCallback(tess, GLU_TESS_END, endCallback)
		gluTessBeginPolygon(tess, 0)

		gluTessBeginContour(tess)
		
		# The coordinates are approximated to integers, since that's how OpenGL understands them.
		for i in self.vertexes:
			if np.isnan(i.x) or np.isnan(i.y): return # triangles aren't updated if there is a NaN value
			point3d = (int(i.x),int(i.y),0)
			gluTessVertex(tess,point3d, point3d) # we put each vertex here
		gluTessEndContour(tess)

		gluTessEndPolygon(tess)
		gluDeleteTess(tess)
		self.triangles = []
		for i in range(0,len(vertexes),3):
			self.triangles.append(vertexes[i:i+3]) # here, we store the triangles' data in the polygon

	# Checks if the vertex defined by x and y can be added to the polygon.
	def canAddVertex(self, x, y):
		# Vertexes cannot be made too close to each other.
		if Vertex(x,y) == self.vertexes[-1]:
			return False
		if len(self.vertexes) > 2:
			for i in range(0,len(self.vertexes)):
				if i+1 >= len(self.vertexes): break
				# If a new side would intersect with another from the same polygon, the new vertex is not added.
				if checkForIntersection(self.vertexes[i], self.vertexes[i+1], self.vertexes[-1], Vertex(x,y)):
					return False
					
		return True
			
	# As the user clicks, more vertexes are added to the polygon. 
	def addVertex(self, x, y):
		if self.canAddVertex(x, y):
			self.vertexes.append(Vertex(x, y))
		
	# Checks if the side made from the first and last vertexes can complete a polygon.
	# A complete polygon needs at least 3 vertexes and cannot have sides intersecting each other.
	def canFinishPolygon(self):
		if len(self.vertexes) > 2:
			for i in range(0,len(self.vertexes)):
				if i+1 >= len(self.vertexes): break
				if checkForIntersection(self.vertexes[i], self.vertexes[i+1], self.vertexes[-1], self.vertexes[0]):
					return False
		else: return False
		
		return True
		
	# Method that moves the polygon to the mouse's position.
	# Each polygon node also moves its children and nails.
	def movePolygon(self, x, y):
		global grabbedPoint
		global grabbedPolygon
		
		# Matrix that translates the polygon. It is shifted in order to move the polygon in relation to
		# the point in which it was first clicked.
		translationMatrix = matrix.translate(-grabbedPoint.x, -grabbedPoint.y, 0) * matrix.translate (x, y, 0)
		
		# The polygon also moves the nail that attached it to a father, if there is one.
		if self.fatherNail != None:
			nailVector = [[self.fatherNail.x], [self.fatherNail.y], [0], [1]]
			translatedNail = np.array(translationMatrix * nailVector)
			if np.isnan(translatedNail[0][0]) or np.isnan(translatedNail[1][0]): return # nail position isn't updated if there is a NaN value
			self.fatherNail.x = translatedNail[0][0]
			self.fatherNail.y = translatedNail[1][0]
			
		# Each vertex of the polygon is translated.
		for vertex in self.vertexes:
			vertexVector = [[vertex.x], [vertex.y], [vertex.z], [1]]
			translatedVector = np.array(translationMatrix * vertexVector)
			vertex.x = translatedVector[0][0]
			vertex.y = translatedVector[1][0]
			
		self.triangulate() # the polygon's respective triangles are updated
		
		# Recursively, each descendant of the polygon is translated as well
		for child in self.children:
			child.movePolygon(x,y)
			
		# The root updates the mouse position that is grabbing the polygon, after all 
		# polygons in the tree have moved.
		if self == grabbedPolygon: grabbedPoint = Vertex(x,y)
		
	# If the polygon is not the root of its tree, it shall be rotated around its father nail.							
	def rotatePolygon(self, x, y):
		global grabbedPoint
		global grabbedPolygon
		global rotationPoint
	
		# Returns the angle between three points, in degrees.
		# Useful for checking how much the polygon should rotate.
		# Points a, b and c must be in the form np.array([x, y, z]).
		def angleBetweenPoints(a, b, c):
			ba = a - b
			bc = c - b

			# Avoids a division by zero
			if (np.linalg.norm(ba) * np.linalg.norm(bc)) == 0: return 0
			
			cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
			
			# Avoids getting the arccos of a nonsensical cosine
			if cosine_angle > 1: cosine_angle = 1
			elif cosine_angle < -1: cosine_angle = -1

			angle = np.arccos(cosine_angle)
			
			if np.cross(ba,bc)[2] > 0:
				return np.degrees(angle)
			else:
				return np.degrees(-angle)
	
		if self == grabbedPolygon: rotationPoint = self.fatherNail
	
		# Matrix that rotates the polygon around the father nail.
		rotationMatrix = matrix.rotate(angleBetweenPoints( \
								np.array([grabbedPoint.x, grabbedPoint.y, 0]), \
								np.array([rotationPoint.x, rotationPoint.y, 0]), \
								np.array([x, y, 0])\
								), \
								0, 0, 1)
								
		# The defacto matrix that rotates the polygon. It shifts the polygon so that the father nail
		# is in the origin, then it rotates around the father nail and finally shifts
		# back to the original position.
		nailRotationMatrix = matrix.translate(rotationPoint.x, rotationPoint.y, 0) * \
				rotationMatrix * \
				matrix.translate(-rotationPoint.x, -rotationPoint.y, 0)
				
		# The polygon also moves the nail that attached it to a father. Useful for the recursion.
		nailVector = [[self.fatherNail.x], [self.fatherNail.y], [0], [1]]
		rotatedNail = np.array(nailRotationMatrix * nailVector)
		if np.isnan(rotatedNail[0][0]) or np.isnan(rotatedNail[1][0]): return # nail position isn't updated if there is a NaN value
		self.fatherNail.x = rotatedNail[0][0]
		self.fatherNail.y = rotatedNail[1][0]
		
		# Each vertex of the polygon is rotated.
		for vertex in self.vertexes:
			vertexVector = [[vertex.x], [vertex.y], [vertex.z], [1]]
			rotatedVector = np.array(nailRotationMatrix * vertexVector)
			vertex.x = rotatedVector[0][0]
			vertex.y = rotatedVector[1][0]
			
		self.triangulate() # the polygon's respective triangles are updated
		
		# Recursively, each descendant of the polygon is translated as well
		for child in self.children:
			child.rotatePolygon(x,y)
			
		# The root updates the mouse position that is grabbing the polygon, after all 
		# polygons in the tree have moved.
		if self == grabbedPolygon: grabbedPoint = Vertex(x,y)

	# Deletes the polygon and removes all nails associated with it.
	def deletePolygon(self):
		global grabbedPolygon
		global grabbedPoint
		
		grabbedPolygon = None
		grabbedPoint = None
	
		if self.father != None:
			self.father.children.remove(self)
			nails.remove(self.fatherNail)
		
		if len(self.children) > 0:
			for child in self.children:
				child.father = None
				nails.remove(child.fatherNail)
				child.fatherNail = None
		
		polygons.remove(self)
		
class Nail(object):
	# A nail is determined by its coordinates, as well as what polygons are being attached by it
	def __init__(self, x, y, father, child):
		self.x = x
		self.y = y
		self.father = father
		self.child = child
		
# Function that manages all commands from the keyboard input.
def keyboardFunc(key, x, y):
	global alpha
	global stroke
	global strokeSize
	global grayScale
	global showNails
	global showTriangles
	global helpActive
	global grabbedPolygon
	global isDrawingPolygon
	global nails
	global polygons

	if key == 't' or key == 'T':
		if not helpActive:
			if alpha == 1: alpha = .75
			else: alpha = 1

	elif key == 'h' or key == 'H':
		helpActive = not helpActive

	# If the user isn't doing anything, C clears the program.
	elif key == 'c' or key == 'C':
		if not helpActive:
			if grabbedPolygon == None and isDrawingPolygon == False:
				del nails[:]
				del polygons[:]
		
	elif key == 's' or key == 'S':
		if not helpActive:
			stroke = not stroke
			
	elif key == '-':
		if not helpActive and stroke:
			if strokeSize > 1: strokeSize -= 1
			
	elif key == '+':
		if not helpActive and stroke:
			if strokeSize < 5: strokeSize += 1

	elif key == 'n' or key == 'N':
		if not helpActive:
			showNails = not showNails

	elif key == 'g' or key == 'G':
		if not helpActive:
			grayScale = not grayScale
			
	elif key == 'i' or key == 'i':
		if not helpActive:
			showTriangles = not showTriangles
			
	elif key == 'r' or key == 'R':
		if not helpActive:
			if grabbedPolygon != None:
				grabbedPolygon.randomizeColor()
			else:
				for p in polygons:
					p.randomizeColor()
					
	elif key == 'd' or key == 'D':
		if not helpActive:
			if grabbedPolygon != None:
				grabbedPolygon.deletePolygon()
		
# Mouse function that adds polygons to the global list.
# Each vertex is defined by the position where the left mouse button is pressed.
# Nails are placed in the position where the right mouse button is pressed,
# as long as there are 2 polygons in that location.
def mouseClick(button, state, x, y):
	global helpActive
	global polygons
	global nails
	global isDrawingPolygon
	global windowHeight
	global mouseX
	global mouseY
	global grabbedPoint
	global grabbedPolygon
	
	# The user can only draw while the help menu isn't shown.
	if not helpActive:
		if (button == GLUT_LEFT_BUTTON and state == GLUT_DOWN):
			if not isDrawingPolygon:
				# Grabs the polygon clicked by the mouse that is closest to the viewer.
				for polygon in reversed(polygons):
					if geometry.Polygon(polygon.vertexes).contains(geometry.Point(x,windowHeight-y)):
						grabbedPolygon = polygon
						grabbedPoint = Vertex(x, windowHeight-y)
						return
						
				# If no polygon was clicked, begins drawing a new polygon.
				polygons.append(Polygon(x, windowHeight-y))
				isDrawingPolygon = True
				
			# If the user is currently drawing a polygon and clicks once more.
			else:
				# Tests whether the current polygon can be finished, and if it can, it becomes a complete polygon and the drawing ends.
				if polygons[-1].vertexes[0] == Vertex(x, windowHeight-y) and len(polygons[-1].vertexes) > 2:
					if polygons[-1].canFinishPolygon():
						isDrawingPolygon = False
						polygons[-1].triangulate()
				# If not, the position the user clicked becomes a new vertex of the polygon, as long as it does not intersect with one of its sides.
				else: polygons[-1].addVertex(x, windowHeight-y)
				
		# When the mouse button is raised, the polygon isn't moved anymore.
		if (button == GLUT_LEFT_BUTTON and state == GLUT_UP):
			if grabbedPolygon != None:
				grabbedPoint = None
				grabbedPolygon = None
				rotationPoint = None
				
		if (button == GLUT_RIGHT_BUTTON and state == GLUT_DOWN):
			if not isDrawingPolygon:
				# Remove a nail in the location
				for nail in nails:
					if Vertex(nail.x, nail.y) == Vertex(x, windowHeight-y):
						nail.father.children.remove(nail.child)
						nail.child.father = None
						nail.child.fatherNail = None
						nails.remove(nail)
						return
						
				# Add a nail to a pair of polygons
				for i in reversed(range(0, len(polygons))):
					if geometry.Polygon(polygons[i].vertexes).contains(geometry.Point(x,windowHeight-y)):
						for j in reversed(range(0, i)):
							if geometry.Polygon(polygons[j].vertexes).contains(geometry.Point(x,windowHeight-y)):
								if polygons[i].father == None:
									polygons[i].father = polygons[j]
									polygons[i].fatherNail = Nail(x, windowHeight-y, polygons[j], polygons[i])
									polygons[j].children.append(polygons[i])
									nails.append(polygons[i].fatherNail)
									return
				
			# Cancel drawing a polygon.
			else:
				isDrawingPolygon = False
				polygons.pop()
							
	glutPostRedisplay()
		

# Mouse function that registers the position of the cursor
# for drawing temporary lines
def mouseMovement(x, y):
	global isDrawingPolygon
	global windowHeight
	global windowWidth
	global mouseX
	global mouseY
	global grabbedPolygon
	
	#if (isDrawingPolygon):
	mouseX = x
	mouseY = windowHeight - y
	
	if grabbedPolygon != None:
		# If grabbing a root polygon
		if grabbedPolygon.father == None: 
			grabbedPolygon.movePolygon(mouseX, mouseY)
		
		# If grabbing a non-root polygon
		else:
			grabbedPolygon.rotatePolygon(mouseX, mouseY)
		
def resizeWindow(w, h):
	global windowHeight
	global windowWidth
	
	# Prevents window from having zero height or width
	if(h == 0):
		h = 1
	if(w == 0):
		w = 1
	
	glMatrixMode(GL_PROJECTION) # sets the current matrix to projection
	glLoadIdentity() # multiply the current matrix by identity matrix
	gluOrtho2D(0, w, 0, h) # sets the camera's size as the full window area
	
	windowHeight = h
	windowWidth = w
	glViewport(0, 0, w, h)
	
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	gluLookAt(0,0,1,0,0,0,0,1,0)
	
# Simple intersection function between two line segments.
def checkForIntersection( p0, p1, p2, p3 ) :
	# Ignores intersections of line segments that start at the same point.
	if p0 == p2 or p0 == p3 or p1 == p2 or p1 == p3:
		return False

	s10_x = p1.x - p0.x
	s10_y = p1.y - p0.y
	s32_x = p3.x - p2.x
	s32_y = p3.y - p2.y

	denom = s10_x * s32_y - s32_x * s10_y

	# Collinear
	if denom == 0:
		return True

	s02_x = p0.x - p2.x
	s02_y = p0.y - p2.y

	s_numer = s10_x * s02_y - s10_y * s02_x

	# No intersection
	if (s_numer < 0) == (denom > 0):
		return False

	t_numer = s32_x * s02_y - s32_y * s02_x

	# No intersection
	if (t_numer < 0) == (denom > 0):
		return False

	# No intersection
	if (s_numer > denom) == (denom > 0) or (t_numer > denom) == (denom > 0):
		return False

	# Intersection detected

	t = t_numer / denom

	intersection_point = [ p0.x + (t * s10_x), p0.y + (t * s10_y) ]

	return intersection_point
	
# Draws a polygon by drawing each of its triangles.
# Necessary for drawing concave polygons, as GL_POLYGON can only draw convex ones.
# Iterates through all of the polygon's triangles and draws them using the polygon's color.
def drawPolygon(polygon):
	global alpha
	global stroke
	global strokeSize
	global grayScale
	global showTriangles
	global showNails

	# Simple function that translates a color to grayscale using the luminosity method
	def toGrayScale(color):
		sum = 0.21 * color[0] + 0.72 * color[1] + 0.07 * color[2]

		return [sum, sum, sum]

	for triangle in polygon.triangles:
		if not grayScale:
			glColor4f(polygon.color[0],polygon.color[1],polygon.color[2], alpha) # sets the color for the polygon as the one stored in it
		else:
			grayColor = toGrayScale(polygon.color)
			glColor4f(grayColor[0],grayColor[1],grayColor[2], alpha) # sets the color for the polygon as the grayscale version of the one stored in it
		glBegin(GL_TRIANGLES)
		for vertex in triangle:
			glVertex2i(int(vertex[0]),int(vertex[1]))
		glEnd()
		
	# Explicitly shows the triangulation of the polygons
	if showTriangles:
		glLineWidth(1) # sets the width of the triangle's internal lines
		for triangle in polygon.triangles:
			if not grayScale:
				glColor4f(polygon.color[0] * .8,polygon.color[1] * .8,polygon.color[2] * .8, alpha) # sets the color for the lines as a slightly darker version of the polygon's
			else:
				grayColor = toGrayScale(polygon.color)
				glColor4f(grayColor[0] * .8,grayColor[1] * .8,grayColor[2] * .8, alpha) # sets the color for the lines as a slightly darker version of the polygon's grayscale
			
			for i in range(0, len(triangle)):
				glBegin(GL_LINES)
				glVertex2i(int(triangle[i][0]), int(triangle[i][1]))
				if i < len(triangle) - 1:
					glVertex2i(int(triangle[i+1][0]), int(triangle[i+1][1]))
				else: # for the very last side of the triangle
					glVertex2i(int(triangle[0][0]), int(triangle[0][1]))
				glEnd()
			

	# Draw the stroke
	if stroke:
		glColor4f(0, 0, 0, alpha) # sets the color for the stroke as black
		glLineWidth(strokeSize) # sets the width of the stroke
		for i in range(0, len(polygon.vertexes)):
			glBegin(GL_LINES)
			glVertex2i(int(polygon.vertexes[i].x), int(polygon.vertexes[i].y))
			if i < len(polygon.vertexes) - 1:
				glVertex2i(int(polygon.vertexes[i+1].x), int(polygon.vertexes[i+1].y))
			else: # for the very last side of the polygon
				glVertex2i(int(polygon.vertexes[0].x), int(polygon.vertexes[0].y))
			glEnd()

	# Draw the nail that is on top of this specific polygon
	if showNails:
		if polygon.fatherNail != None:
			glColor3f(.3,.3,.3) # sets the color for the nails as grey
			glBegin(GL_POINTS)
			glVertex2i(int(polygon.fatherNail.x), int(polygon.fatherNail.y)) # the coordinates are approximated to integers, since that's how OpenGL understands them
			glEnd()
	
# Function that draws the help text.
def drawHelp(text, x, y):
	# Increases the spacing between each line.
	line = 0
	# Increases the spacing in the beginning of the line.
	indent = 0

	glPushMatrix()
	glLoadIdentity()
	glRasterPos2i(x,y)
	for i in range(0, len(text)):
		if text[i] == str(""): continue # ignores the \ character used to divide lines in Python
		if text[i] == str("\n"):
			line += 20
			indent = 0
			glRasterPos2i(x,y - line)
		elif text[i] == str("\t"):
			indent += 20
			glRasterPos2i(x + indent,y - line)
		else:
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(text[i]))
	glPopMatrix()
	
# Returns a string with indicators for all the options currently active.
def activeOptions():
	global alpha
	global stroke
	global strokeSize
	global grayScale
	global showNails
	global showTriangles
	
	string = ""
	
	if grayScale: string = string + "G"
	else: string = string + " "
	
	if alpha < 1: string = string + "T"
	else: string = string + " "
	
	if showNails: string = string + "N"
	else: string = string + " "
	
	if showTriangles: string = string + "I"
	else: string = string + " "
	
	if stroke: string = string + "S" + "(" + str(strokeSize) + ")"
	else: string = string + " "
	
	return string

# Function that draws the scene every frame.
def renderScene():
	global windowHeight
	global windowWidth
	global polygons
	global nails
	global lineIntersections
	global mouseX
	global mouseY
	global helpActive
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); # clears the buffer to draw new frame

	# Shows the help menu, and only that.
	if helpActive:
		# Draws a vertical gradient background for the help menu
		glBegin(GL_QUADS)
		# White on top
		glColor3f(1.0,1.0,1.0)
		glVertex3f(0.0, windowHeight, 0)
		glVertex3f(windowWidth,windowHeight, 0)
		# Light blue on bottom
		glColor3f(0.8,0.8,1.0)
		glVertex3f(windowWidth,0.0, 0)
		glVertex3f(0.0, 0.0, 0)
		glEnd()
		
		glColor3f(0,0,0) # sets the color for the help indicator 
		drawHelp("Author: Julio Rama \
				\n \
				\nMade with: \
				\n\tPython 2.7 \
				\n\tPyOpenGL \
				\n\tNumpy \
				\n\tgeometry.py \
				\n\tmatrix.py \
				\n \
				\nCommands: \
				\n\tH: Toggle help menu. \
				\n\tG: Toggle grayscale. \
				\n\tT: Toggle transparency. \
				\n\tN: Hide/Show nails. \
				\n\tI: Toggle explicit triangulation of polygons. \
				\n\tS: Toggle stroke. \
				\n\t\tIf stroke is active: \
				\n\t\t+: Increase stroke size. \
				\n\t\t-: Decrease stroke size. \
				\n\tC: Clear screen. \
				\n\tD: Delete a polygon being held. \
				\n\tR: Randomize color of a polygon being held, or all polygons if none is being held. \
				\n \
				\n\tLeft Mouse Button: Draw or hold polygons. Dragging moves the held polygon. \
				\n\tRight Mouse Button: Place nails on polygons or cancel drawing one. \
				\n \
				\nPolygons cannot have intersecting sides.",\
			20, windowHeight - 20)
	
	else:
		if len(polygons) > 0:
			for i in range(0,len(polygons)):
				if isDrawingPolygon and i == len(polygons) - 1:
					break
				drawPolygon(polygons[i])
		
		glColor3f(0,0,0) # sets the color for the lines as black
		glLineWidth(1) # sets the width of the lines of the new polygon
		if isDrawingPolygon:
			lastPolygon = polygons[-1]
			for i in range(0,len(lastPolygon.vertexes)):
				if i+1 >= len(lastPolygon.vertexes): break
				glBegin(GL_LINES)
				glVertex2i(lastPolygon.vertexes[i].x,lastPolygon.vertexes[i].y)
				glVertex2i(lastPolygon.vertexes[i+1].x,lastPolygon.vertexes[i+1].y)
				glEnd()
				
			if not lastPolygon.canAddVertex(mouseX, mouseY):
				glColor4f(.5,.5,.5, .5)
						
			glBegin(GL_LINES)
			glVertex2i(lastPolygon.vertexes[-1].x,lastPolygon.vertexes[-1].y)
			glVertex2i(mouseX,mouseY)
			glEnd()
			
			if lastPolygon.canFinishPolygon() and len(lastPolygon.vertexes) > 2 and lastPolygon.vertexes[0] == Vertex(mouseX,mouseY):
				glColor4f(0,0,0,.5) # sets the color for the point as transparent black 
				glBegin(GL_POINTS)
				glVertex2i(lastPolygon.vertexes[0].x, lastPolygon.vertexes[0].y)
				glEnd()
				
		glColor3f(0.4,0.4,0.4) # sets the color for the help indicator 
		drawHelp("H: help menu", \
				20, windowHeight - 20)
				
		
		glColor3f(0.4,0.4,1.0) # sets the color for the active options indicator 
		drawHelp("Active: " + activeOptions(), \
				20, 20)
		
	glFlush()
	

def main():
	glutInit(sys.argv)
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)
	glutInitWindowPosition(100,100)
	glutInitWindowSize(800,600)
	glutCreateWindow("Trabalho 2 CG - Julio Rama")
	glClearColor(1, 1, 1, 1) # sets the background color to white
	glEnable(GL_POINT_SMOOTH) # makes the points circular instead of squares
	glutDisplayFunc(renderScene)
	glPointSize(10) # sets the size of the points and nails
	glEnable(GL_BLEND) # enables blending for transparent lines
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) # sets the blending function
	glutIdleFunc(renderScene)
	glutKeyboardFunc(keyboardFunc)
	glutMouseFunc(mouseClick)
	glutMotionFunc(mouseMovement)
	glutPassiveMotionFunc(mouseMovement)
	glutReshapeFunc(resizeWindow)
	glutMainLoop()

if __name__ == '__main__': main()
