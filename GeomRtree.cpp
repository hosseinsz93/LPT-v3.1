#include "GeomRtree.h"
#include "GeomSpatialIndexSource.h"

#include "GeomSpatialIndexQuery.h"
#include "GeomDistanceSpatialIndexQuery.h"
#include "GeomRaySpatialIndexQuery.h"

#include "GeomException.h"

#include "GeomBoundingBox.h"
#include "GeomCheckInterrupts.h"

#include "GeomId.h"
#include "GeomSubset.h"

#include "GeomBox3.h"
#include "GeomPoint3.h"

#include <limits>
#include <algorithm> // for sorting
#include <queue> // for priority queue
#include <vector>

namespace Geom
{
// smallest number greater than one
static float const ONE_PLUS_EPSILON =
    1.0f + std::numeric_limits<float>::epsilon();

// largest number less then one
static float const ONE_MINUS_EPSILON =
    1.0f - std::numeric_limits<float>::epsilon();

/**
 * Convert a double to a float, nudging it towards zero.  The goal here is
 * to guarantee that on conversion back to a double (which happens during
 * queries), the new double evaluates <= the original double.
 */
inline float convertDownToFloat(double d) {
  float f = d;
  return (f > 0 ? f*ONE_MINUS_EPSILON : f*ONE_PLUS_EPSILON);
}

/**
 * Convert a double to a float, nudging it away from zero.  The goal here is
 * to guarantee that on conversion back to a double (which happens during
 * queries), the new double evaluates >= the original double.
 */
inline float convertUpToFloat(double d) {
  float f = d;
  return (f > 0 ? f*ONE_PLUS_EPSILON : f*ONE_MINUS_EPSILON);
}

/*! An internal class used to store an MBC (Minimal Bounding Cube)
    of an object. */
struct RtreeBox
{
  /*! Default constructor -- min is set to huge numbers and max is
      set to huge negative numbers so the first call to expand will
      initialize box. */
  RtreeBox() {
    this->clear();
  }

  inline void clear() {
    min(0) = std::numeric_limits<float>::max();
    min(1) = min(0);
    min(2) = min(0);

    max(0) = -min(0);
    max(1) = max(0);
    max(2) = max(0);
  }

  /*! Construct from a BoundingBox object. */
  RtreeBox(BoundingBox const &box) {
    // inflate bounding box slightly on creation.  Since we're converting
    // double precision coordinates down to floats then possibly back up
    // to doubles for queries, we need to nudge the bounding box coordinates
    // a bit to make up for loss of precision.
    min(0) = convertDownToFloat(box.min()(0));
    min(1) = convertDownToFloat(box.min()(1));
    min(2) = convertDownToFloat(box.min()(2));

    max(0) = convertUpToFloat(box.max()(0));
    max(1) = convertUpToFloat(box.max()(1));
    max(2) = convertUpToFloat(box.max()(2));
  }

  /*! Expand to include \a box. */
  void expand(RtreeBox const &box) {
    min(0) = std::min(min(0), box.min(0));
    min(1) = std::min(min(1), box.min(1));
    min(2) = std::min(min(2), box.min(2));

    max(0) = std::max(max(0), box.max(0));
    max(1) = std::max(max(1), box.max(1));
    max(2) = std::max(max(2), box.max(2));
  }

  /*! Compute and return the center of the box. */
  void getCenter(float *center) const {
    center[0] = 0.5*(min(0) + max(0));
    center[1] = 0.5*(min(1) + max(1));
    center[2] = 0.5*(min(2) + max(2));
  }

  /*! Compute and return the longest dimension (0, 1, or 2) of the
      box. */
  int getLongestDimension() const {
    float dx = max(0) - min(0);
    float dy = max(1) - min(1);
    float dz = max(2) - min(2);

    if (dx > dy)
    {
      if (dx > dz)
        return 0;
      return 2;
    }
    if (dy > dz)
      return 1;
    return 2;
  }

  /*! Returns true if intersects \a box, false otherwise */
  bool intersects(RtreeBox const &box) const {
    return (min(0) <= box.max(0) && max(0) >= box.min(0) &&
        min(1) <= box.max(1) && max(1) >= box.min(1) &&
        min(2) <= box.max(2) && max(2) >= box.min(2));
  }


  /*! Returns true if intersects \a box, false otherwise */
  bool intersects(Box3 const &box) const {
    return (min(0) <= box.max()(0) && max(0) >= box.min()(0) &&
        min(1) <= box.max()(1) && max(1) >= box.min()(1) &&
        min(2) <= box.max()(2) && max(2) >= box.min()(2));
    return true;
  }

  /*! Returns true if the box is valid (has been resized from it's initial
      size), false otherwise. */
  bool isValid() const {
    return (max(0) > min(0));
  }

  /*! Compute the volume of the box. */
  float volume() const {
    if (!this->isValid())
      // box must be empty -- return zero volume
      return 0.0;

    return ((max(0) - min(0)) *
        (max(1) - min(1)) *
        (max(2) - min(2)));
  }

  float coverage() const {
    return volume();
  }

  float computeMargin() const {
    return max(0)-min(0) +max(1)-min(1) + max(2)-min(2);
  }
  /*! added solely for purpose of allowing RtreeBox to become
      a virtual class so that 2d boxes and 3d boxes could be used
      and the insertion/deletion algorithm can still work
      just returns the dimensionality of this box
   */
  unsigned getNumDimensions() const{
    return 3;
  }
  /*! returns min[axis] either the x,y or z value of the min coord
   */
  float getMinAxisValue(int axis) const{
    return min(axis);
  }
  /*! returns max[axis] either the x,y or z value of the max coord
   */
  float getMaxAxisValue(int axis) const{
    return max(axis);
  }
  /*!
      if this box intersects with box2 this method
      computes intersection of this box with box2 and puts the overlap
      in intsectBox and returns true- 
      if no intersection returns false and intsectBox is undefined
   */
  bool computeIntersection(RtreeBox const &box2,RtreeBox &intsectBox)const{
    for (unsigned i=0;i<box2.getNumDimensions();i++)
    {
      intsectBox.min(i)=std::max(this->min(i),box2.min(i));
      intsectBox.max(i)=std::min(this->max(i),box2.max(i));
      if (intsectBox.min(i)>=intsectBox.max(i))
        return false;
    }
    return true;
  }
  /*!   returns true if box2 is completely contained within this box
      false otherwise
   */
  bool contains(RtreeBox const &box2) const {
    return ((box2.min(0) >= this->min(0))&&
        (box2.min(1) >= this->min(1))&&
        (box2.min(2) >= this->min(2))&&
        (box2.max(0) <= this->max(0))&&
        (box2.max(1) <= this->max(1))&&
        (box2.max(2) <= this->max(2)));
  }
  bool operator==(RtreeBox const &rhs)const
                    {
    return (rhs.min(0)==min(0) &&
        rhs.min(1)==min(1) &&
        rhs.min(2)==min(2) &&
        rhs.max(0)==max(0) &&
        rhs.max(1)==max(1) &&
        rhs.max(2)==max(2));
                    }
  friend std::ostream &operator<<(std::ostream &loutstr,
                                  RtreeBox const &box)
  {
    loutstr << "("
            << box.min(0) << ","
            << box.min(1) << ","
            << box.min(2)
            << ")-"
            << "("
            << box.max(0) << ","
            << box.max(1) << ","
            << box.max(2)
            << ")";
    return loutstr;
  }

    Vector3f min;
    Vector3f max;
};


// forward declare needed for search queue
class RtreeNode;

/* Internal Stuff for priority queue for searching */
class RtreeSearchQueueObject
{
public:
  RtreeSearchQueueObject()
  : _node(0)
  , _distance(0.0)
  {}

  RtreeSearchQueueObject(RtreeNode const *node, double distance)
  : _node(node)
  , _distance(distance)
  {}

  RtreeSearchQueueObject(Id const &id, double distance)
  : _node(0)
  , _id(id)
  , _distance(distance)
  {}

  bool isNode() const {return _node != 0;}
  bool isElement() const {return _id.isValid();}
  RtreeNode const *getNode() const {return _node;}
  Id const &getId() const {return _id;}
  double getDistance() const {return _distance;}

  bool operator<(RtreeSearchQueueObject const &rhs) const {
    //fixme? make this check with a tolerance?
    if (_distance > rhs._distance)
      return true;
    else if (_distance < rhs._distance)
      return false;

    // distances are the same -- check type.  Need to make sure that
    // id objects get sorted to the left of node objects with the same
    // distance.  Sorting between objects of the same type is irrelevant,
    // so always return false in this case.

    if (this->isNode())
      return false; // lhs is a node.  sort to the right regardless of
    // what RHS is

    if (rhs.isNode())
      return true; // lhs is an element, rhs is a node.  Sort lhs to the left

    //both are elements, order by id so ordering is consistent
    return rhs.getId() < this->getId();
  }
private:
  RtreeNode const *_node;
  Id _id;
  double _distance;
};


class RtreeSearchQueue
{
public:
  RtreeSearchQueue() {}
  ~RtreeSearchQueue() {}

  // clear the queue
  void clear() {
    _queue = QueueType();
  }

  // check if the queue is empty
  bool empty() const {return _queue.empty();}

  // push an Id onto the queue
  void pushId(Id const &id, double distance) {
    _queue.push(RtreeSearchQueueObject(id, distance));
  }
  // push a node onto the queue
  void pushNode(RtreeNode const *node, double distance) {
    _queue.push(RtreeSearchQueueObject(node, distance));
  }

  // pop the last object off the queue
  void pop() {
    _queue.pop();
  }

  // retrieve the top object off the queue.
  RtreeSearchQueueObject top() const {
    return _queue.top();
  }
private:
  typedef std::priority_queue<RtreeSearchQueueObject> QueueType;
  QueueType _queue;
};


/*! An internal class to store an Id with its bounding box in an R-tree. */
struct RtreeLeafEntry
{
  RtreeLeafEntry(RtreeBox const &box_, Id const &id_)
  : box(box_), id(id_) {}
  RtreeLeafEntry(BoundingBox const &box_, Id const &id_)
  : box(box_), id(id_) {}
  RtreeLeafEntry(){};
  inline void recomputeBoundingBox(SpatialIndexSource const &/*source*/) {

  }

  RtreeBox box;
  Id id;
};

/*! A base class storing an R-tree node, could be leaf or non-leaf. */
struct RtreeNode
{
  RtreeNode() {}
  virtual ~RtreeNode() {}

  /*! Add all objects who's MBR intersects \a box to vector \a obj. */
  virtual void addIntersectingObjects(RtreeBox const &box,
                                      std::vector<Id> &obj,
                                      SpatialIndexSource const &source) = 0;

  /*! Check for any objects who's MBR intersects \a box  */
  virtual bool checkIntersectingObjects(RtreeBox const &box,
                                        SpatialIndexSource const &source) = 0;

  /*! Enqueue the children of this node into a search queue with a given
      query function. */
  virtual void enqueueChildren(RtreeSearchQueue &searchqueue,
                               SpatialIndexQuery *query) const = 0;

  virtual void recomputeBoundingBox(SpatialIndexSource const &source) = 0;

  virtual bool recomputeBoundingBoxContaining(Id const &id,
                                              SpatialIndexSource const &source) = 0;

  virtual bool recomputeBoundingBoxContaining(Subset const &subset,
                                              SpatialIndexSource const &source) = 0;
  virtual bool isLeaf()const = 0;

  RtreeBox const & getBox()const {return box;};
  virtual RtreeBox const & getChildBox(unsigned index) const = 0;
  virtual unsigned getNumChildren() const = 0;
  virtual void removeChild(unsigned childIndex) = 0;
  virtual void clearChildren()=0;
  virtual void print(bool decend = true) const=0;
  virtual void print(int level, bool decend = true) const=0;
  virtual void printBox(int level) const
  {
    for (int i=0; i<level; ++i)
      std::cout << "  ";
    std::cout << box << std::endl;
  }
  RtreeBox box;
};

/*! Non-leaf R-tree node. */
struct RtreeNonLeafNode : public RtreeNode
{
  /*! Construct the node.  The children object will be reserved to size
      \a numchildren. */
  RtreeNonLeafNode(int numchildren) {
    _children.reserve(numchildren);
  }
  RtreeNonLeafNode() {
  }

  ~RtreeNonLeafNode() {
    for (unsigned i = 0; i<_children.size(); ++i)
      delete _children[i];
  }

  RtreeNode *getChild(unsigned index) {
    return _children[index];
  }
  void addIntersectingObjects(RtreeBox const &bbox,
                              std::vector<Id> &obj,
                              SpatialIndexSource const &source)
  {
    // recurse to children
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i]->box.intersects(bbox))
        _children[i]->addIntersectingObjects(bbox, obj, source);
    }
  }

  bool checkIntersectingObjects(RtreeBox const &bbox,
                                SpatialIndexSource const &source)
  {
    // recurse to leaf children
    for(unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i]->box.intersects(bbox))
        if(_children[i]->checkIntersectingObjects(bbox, source))
          return true;
    }
    return false;
  }

  /*! Add a child node.  The child nodes MBR must already be computed. */
  void addChild(RtreeNode *node) {
    _children.push_back(node);
    box.expand(node->box);
  }

  void enqueueChildren(RtreeSearchQueue &searchqueue,
                       SpatialIndexQuery *query) const {
    static Box3 scratchBox;

    double distance;
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      RtreeNode *child = _children[i];

      if (query->hasQueryBox())
      {
        // do (faster) overlap test before (more expensive) distance test
        if (!child->box.intersects(query->getQueryBox()))
          // overlap test failed
          continue;
      }
      scratchBox.min() = child->box.min;
      scratchBox.max() = child->box.max;
      if (query->getDistance(scratchBox,
          &distance))
        searchqueue.pushNode(child, distance);
    }
  }

  void recomputeBoundingBox(SpatialIndexSource const &source) {
    box.clear();
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      _children[i]->recomputeBoundingBox(source);
      box.expand(_children[i]->box);
    }
  }

  bool recomputeBoundingBoxContaining(Id const &id,
                                      SpatialIndexSource const &source) {
    bool recomputeMyBox = false;
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i]->recomputeBoundingBoxContaining(id, source))
        recomputeMyBox = true;
    }

    if (recomputeMyBox)
    {
      box.clear();
      for (unsigned i = 0; i<_children.size(); ++i)
        box.expand(_children[i]->box);
      return true;
    }
    return false;
  }

  bool recomputeBoundingBoxContaining(Subset const &subset,
                                      SpatialIndexSource const &source) {
    bool recomputeMyBox = false;
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i]->recomputeBoundingBoxContaining(subset,
          source))
        recomputeMyBox = true;
    }

    if (recomputeMyBox)
    {
      box.clear();
      for (unsigned i = 0; i<_children.size(); ++i)
        box.expand(_children[i]->box);
      return true;
    }
    return false;
  }

  RtreeBox const &getChildBox(unsigned index) const {
    return _children[index]->box;
  }
  bool isLeaf() const {
    return false;};
  unsigned getNumChildren() const {
    return _children.size();};
  void clearChildren() {
    _children.clear();
    box.clear();
  }
  void removeChild(unsigned childIndex) {
    _children.erase(_children.begin()+childIndex);
    box.clear();
    for (unsigned i = 0; i<_children.size(); i++)
    {
      box.expand(_children[i]->box);
    }
  }
  bool findChild(RtreeNode *child,unsigned &index) const {
    for (unsigned i=0;i<_children.size();i++)
    {
      if (_children[i]==child)
      {
        index = i;
        return true;
      }
    }
    return false;
  }
  void updateBox() {
    this->box.clear();
    for (unsigned i=0;i<_children.size();i++)
      this->box.expand(_children[i]->box);
  }

  void print(bool decend) const { print(0, decend); }

  void print(int level, bool decend) const
  {
    for (int i=0; i<level; ++i)
      std::cout << "  ";
    std::cout << "RtreeNonLeafNode: " << this << std::endl;
    printBox(level);

    if (decend)
      for (int i=0; i<(int)_children.size(); ++i)
        _children[i]->print(level+1, decend);
  }

private:
  std::vector<RtreeNode *> _children;
};

/*! Leaf R-tree node. */
struct RtreeLeafNode : public RtreeNode
{
  /*! Construct the node.  The children object will be reserved to size
      \a numchildren. */
  RtreeLeafNode(int numchildren) {
    _children.reserve(numchildren);
  }
  RtreeLeafNode() {
  }
  ~RtreeLeafNode() {}

  void addIntersectingObjects(RtreeBox const &bbox,
                              std::vector<Id> &obj,
                              SpatialIndexSource const &)
  {
    // add children to list if their MBR intersects the input box
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i].box.intersects(bbox))
        obj.push_back(_children[i].id);
    }
  }

  bool checkIntersectingObjects(RtreeBox const &bbox,
                                SpatialIndexSource const &)
  {
    // check children for their MBR intersecting the input box
    for(unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i].box.intersects(bbox))
        return true;
      //obj.push_back(_children[i].id);
    }
    return false;
  }



  /*! Add a child entry. */
  void addChild(RtreeLeafEntry const &entry) {
    _children.push_back(entry);
    box.expand(entry.box);
  }

  void enqueueChildren(RtreeSearchQueue &searchqueue,
                       SpatialIndexQuery *query) const {
    double distance;
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      Id child = _children[i].id;
      if (query->getDistance(child,&distance))
        searchqueue.pushId(child, distance);
    }
  }

  void recomputeBoundingBox(SpatialIndexSource const &source) {
    box.clear();
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      _children[i].box = source.getBoundingBox(_children[i].id);
      box.expand(_children[i].box);
    }
  }

  bool recomputeBoundingBoxContaining(Id const &id,
                                      SpatialIndexSource const &source) {
    bool recomputeMyBox = false;
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      if (_children[i].id == id)
      {
        _children[i].box = source.getBoundingBox(id);
        recomputeMyBox = true;
        break; // id can only be contained once
      }
    }
    if (recomputeMyBox)
    {
      box.clear();
      for (unsigned i = 0; i<_children.size(); ++i)
        box.expand(_children[i].box);
      return true;
    }
    return false;
  }

  bool recomputeBoundingBoxContaining(Subset const &subset,
                                      SpatialIndexSource const &source) {
    bool recomputeMyBox = false;
    for (unsigned i = 0; i<_children.size(); ++i)
    {
      if (subset.contains(_children[i].id))
      {
        _children[i].box = source.getBoundingBox(_children[i].id);
        recomputeMyBox = true;
      }
    }
    if (recomputeMyBox)
    {
      box.clear();
      for (unsigned i = 0; i<_children.size(); ++i)
        box.expand(_children[i].box);
      return true;
    }
    return false;
  }

  RtreeLeafEntry const &getChild(unsigned i) const {
    return _children[i];
  }
  void clearChildren() {
    _children.clear();
    box.clear();
  }
  bool isLeaf() const {
    return true;
  }
  RtreeBox const &getChildBox(unsigned index) const {
    return _children[index].box;
  }
  unsigned getNumChildren() const {
    return _children.size();
  }
  bool findEntry(Id entry,unsigned &index) {
    for (unsigned i=0;i<_children.size();i++)
      if (_children[i].id==entry)
      {
        index = i;
        return true;
      }
    return false;
  }
  void removeChild(unsigned childIndex) {
    _children.erase(_children.begin()+childIndex);
    box.clear();
    for (unsigned i = 0; i<_children.size(); ++i)
      box.expand(_children[i].box);
  }

  void print(bool decend) const { print(0, decend); }

  void print(int level, bool decend) const
  {
    for (int i=0; i<level; ++i)
      std::cout << "  ";
    std::cout << "RtreeLeafNode: " << this << std::endl;
    printBox(level);

    for (int i=0; i<level+1; ++i)
      std::cout << "  ";
    for (int i=0; i<(int)_children.size(); ++i)
    {
      if (i != 0)
        std::cout << ",";
      std::cout << _children[i].id;
    }
    std::cout << std::endl;
  }

private:
  std::vector<RtreeLeafEntry> _children;
};

/*! Internal insulation class. */
struct Rtree::Data
{
  RtreeNode *root;
  SpatialIndexQuery *query;
  RtreeSearchQueue searchqueue;
  unsigned minN,maxN,reInsP; //min and max
};

/*! Construct the R-tree based on SpatialIndexSource \a source. */
Rtree::Rtree(std::shared_ptr<SpatialIndexSource> const & source)
: SpatialIndex(source)
{
  _data = new Data();
  _data->root = 0;
  _data->query = 0;
  _data->minN = 3;
  // hard code N for now -- 8 seems like a reasonable number here
  _data->maxN = 8;
  _data->reInsP=unsigned(ceil(0.3*_data->maxN));
}

Rtree::~Rtree()
{
  delete _data->root;
  delete _data;
}

template<class T, class U> T assert_cast(U u) {
  T t = dynamic_cast<T>(u);
  if (!t)
    throw Exception("Dynamic cast error in Rtree");
  return t;
}

// forward declare static method
static RtreeNode *TGSBulkLoad(std::vector<RtreeLeafEntry> const &entries,
                              int N);

void Rtree::addSubset(Subset const &subset)
{
  //! \todo Add 'Small Tree/Large Tree' algorithm to R-tree
  if (_data->root)
    throw Exception("Cannot yet append a second subset onto R-tree");

  CheckInterrupts interrupt(CheckInterrupts::high);
  int const N = _data->maxN;

  // test if input number of elements is less than N.
  // If so, root is a leaf node.
  if ((int)subset.getNumElements() <= N)
  {
    RtreeLeafNode *node = new RtreeLeafNode(subset.getNumElements());
    for (SubsetIterator si(&subset); !si.done(); ++si) {
      interrupt.check();
      node->addChild(RtreeLeafEntry(this->getSource().getBoundingBox(*si),
          (*si)));
    }
    _data->root = node;
    return;
  }

  // allocate a vector of RtreeLeafEntries and fill from subset
  std::vector<RtreeLeafEntry> leafentries;
  leafentries.reserve(subset.getNumElements());

  for (SubsetIterator si(&subset); !si.done(); ++si)
  {
    interrupt.check();
    BoundingBox bbox = this->getSource().getBoundingBox(*si);
    leafentries.push_back(RtreeLeafEntry(bbox, (*si)));
  }

  // call bulk loader.
  _data->root = TGSBulkLoad(leafentries, N);
  interrupt.forceCheck();
}

void Rtree::clear()
{
  delete _data->root;
  _data->root = 0;
}

void Rtree::getOverlappingEntities(BoundingBox const &bbox,
                                   std::vector<Id> &entities) const
{
  RtreeBox box(bbox);
  _data->root->addIntersectingObjects(bbox, entities, this->getSource());
}


bool Rtree::checkOverlappingEntities(BoundingBox const &bbox) const
{
  RtreeBox box(bbox);
  return _data->root->checkIntersectingObjects(bbox, this->getSource());
}


BoundingBox Rtree::getBoundingBox() const
{
  return BoundingBox(Vector3d(_data->root->box.min(0),
      _data->root->box.min(1),
      _data->root->box.min(2)),
                     Vector3d(_data->root->box.max(0),
          _data->root->box.max(1),
          _data->root->box.max(2)));
}

  void Rtree::beginPointQuery(Vector3d const &point,
                              double maxDistance) const
{
  if (_data->query)
    throw Exception("Rtree::beginQuery: Another query still open");

  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");

  this->beginQuery(new DistanceSpatialIndexQuery(this->getSource(),
      point,
      maxDistance));
}

void Rtree::beginRayQuery(Ray3 const &ray,
                          double maxDistance) const
{
  if (_data->query)
    throw Exception("Rtree::beginQuery: Another query still open");

  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");

  this->beginQuery(new RaySpatialIndexQuery(this->getSource(),
      ray.clone(),
      maxDistance));
}

void Rtree::beginQuery(SpatialIndexQuery *query) const
{
  if (_data->query)
    throw Exception("Rtree::beginQuery: Another query still open");

  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");
  _data->query = query;

  // clear priority queue
  _data->searchqueue.clear();

  // push root node onto priority queue
  double distance;
  if (_data->query->getDistance(Box3(Vector3d(_data->root->box.min),
                                     Vector3d(_data->root->box.max)),
                                             &distance))
    _data->searchqueue.pushNode(_data->root, distance);
}


Id Rtree::query(double *object_distance) const
{
  if (!_data->query)
    throw Exception("Rtree::query: No query open");

  while (!_data->searchqueue.empty())
  {
    // retrieve top object from queue and pop it off
    RtreeSearchQueueObject obj = _data->searchqueue.top();
    _data->searchqueue.pop();

    // if (obj.getNode())
    //   obj.getNode()->print(0, false);

    if (obj.isElement())
    {
      // front object is an element -- return it
      if (object_distance)
        *object_distance =
            _data->query->convertDistance(obj.getDistance());
      return obj.getId();
    }

    // object is a node, enqueue its children
    RtreeNode const *node = obj.getNode();
    if (node)
      node->enqueueChildren(_data->searchqueue,
          _data->query);
  }

  return Id();
}

void Rtree::endQuery() const
{
  delete _data->query;
  _data->query = 0;
  _data->searchqueue.clear();
}

  Id Rtree::getClosestEntity(Vector3d const &point,
                             double *object_distance,
                             double maxDistance) const
{
  if (_data->query)
    throw Exception("Rtree::beginQuery: Another query still open");

  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");

  this->beginQuery(new DistanceSpatialIndexQuery(this->getSource(),
      point,
      maxDistance,
      true));
  Id result = this->query(object_distance);
  this->endQuery();
  return result;
}

void Rtree::checkForPrevQuery() const
{
  if (_data->query)
    throw Exception("Rtree::beginQuery: Another query still open");

  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");
}

void Rtree::recomputeAllBoundingBoxes()
{
  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");

  _data->root->recomputeBoundingBox(this->getSource());
}

void Rtree::recomputeBoundingBoxesContaining(Id const &id)
{
  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");

  _data->root->recomputeBoundingBoxContaining(id, this->getSource());
}

void Rtree::recomputeBoundingBoxesContaining(Subset const &subset)
{
  if (!_data->root)
    throw Exception("Rtree::beginQuery: Nothing in tree");

  if (!subset.empty())
    _data->root->recomputeBoundingBoxContaining(subset, this->getSource());
}

static bool insertEntry(RtreeNode *&root,
                        RtreeLeafEntry &newEntry,
                        Rtree::Data *data,
                        SpatialIndexSource const &source);
static bool removeEntry(RtreeNode *&root,
                        RtreeLeafEntry &delEntry,
                        Rtree::Data *data,
                        SpatialIndexSource const &source,
                        bool assumeValid);

static bool removeEntities(RtreeNode *node,
                           Subset const &deletedEntities,
                           unsigned numEntitiesToDelete,
                           unsigned numDeletedEntities);


bool Rtree::addElement(Id item)
{
  if (_data->query)
    throw Exception("Cannot add elements during a rtree query");
  RtreeLeafEntry newEntry(this->getSource().getBoundingBox(item), item);
  return insertEntry(_data->root, newEntry, _data, this->getSource());
}
bool Rtree::removeElement(Id item)
{
  if (_data->query)
    throw Exception("Cannot remove elements during a rtree query");
  RtreeLeafEntry delEntry(this->getSource().getBoundingBox(item), item);
  return removeEntry(_data->root, delEntry, _data, this->getSource(), true);
}

bool Rtree::removeInvalidElement(Id item)
{
  BoundingBox dummyBox(Vector3d(0.0, 0.0, 0.0));
  RtreeLeafEntry delEntry(dummyBox, item);
  return removeEntry(_data->root, delEntry, _data, this->getSource(), false);
}

void Rtree::removeSubset(Subset const &deletedEntities)
{
  unsigned numEntitiesToDelete = deletedEntities.getNumElements();
  if (numEntitiesToDelete == 0)return; // nothing to do

  unsigned numDeletedEntities = 0;
  removeEntities(_data->root, deletedEntities,
                 numEntitiesToDelete, numDeletedEntities);

  if (_data->root->getNumChildren() == 0)
  {
    // special case -- everything deleted!   Turn root into a leaf node
    delete _data->root;
    _data->root = new RtreeLeafNode();
  }

  //     // now adjust tree
  //     std::vector<Id> orphanedEntries;
  //     removeUnderfullNodes(root, data->minN);
}

void Rtree::print(bool decend) const
{
  _data->root->print(0, decend);
}

static void getAllEntriesInternal(RtreeNode *root,Subset &items,
                                  Subset &duplicates);

void Rtree::getAllEntries(Subset &items,
                          Subset &duplicates)
{
  items.clear();
  duplicates.clear();
  getAllEntriesInternal(_data->root,items,duplicates);
}

// this is just a debugging utility
void Rtree::checkVsSubset(Subset const &ids)
{
  Subset items,duplicates;
  getAllEntries(items,duplicates);
  if (duplicates.getNumElements())
    std::cout << "rtree contains duplicates="
              <<duplicates.getNumElements()
              << ", " << duplicates
              <<std::endl;

  std::cout << "rtree contains " << items.getNumElements()
            << " faces and active contains "
            << ids.getNumElements() << std::endl;

  Subset notInTree;
  for (SubsetIterator si(&ids);!si.done();++si)
  {
    if (!items.contains(*si))
      notInTree.add(*si);
  }

  std::cout << "surface contains " << notInTree.getNumElements()
            << " faces not in rtree: "
            << notInTree
            << std::endl;
  items.remove(ids);

  std::cout << "rtree contains " << items.getNumElements() << " elements "
            << "not in surface: "
            << items
            << std::endl;
}

// the following are ifdef'd out but can be useful for debugging
#ifdef InsertionDebugOn
#include "mk/common/LogStream.h"
//chrisgdebug code
struct LevelStats {
public:
  LevelStats():numNodes(0),totNumChildren(0),childoverlap(0.0),
  margin(0.0),coverage(0.0)
  {};
  std::vector<std::pair<RtreeNode*,unsigned> >nodeList;
  unsigned numNodes,totNumChildren;
  double childoverlap,margin,coverage;
};
typedef std::vector<LevelStats> LevelStatList;

static bool checkBoundingBoxesInternal(RtreeNode *root,
                                       RtreeBox &parentBox,
                                       int level,
                                       Rtree::Data *data,
                                       SpatialIndexSource const &source);
static bool checkLeafLevelsInternal(RtreeNode *root);
static bool findNodeLevel(RtreeNode *root, RtreeNode *node,
                          unsigned &level);
static unsigned getLeafLevel(RtreeNode *start,unsigned level);
static double computeMargin(RtreeNode *root,
                            SpatialIndexSource const &source);
static double computeOverlap(RtreeNode *root,
                             SpatialIndexSource const &source);
static void computeLevelStatsInt(RtreeNode *root,LevelStatList &statList,
                                 unsigned level,
                                 SpatialIndexSource const &source);

void Rtree::printLevelStats()
{
  LevelStatList statList;
  computeLevelStatsInt(_data->root,statList,0,this->getSource());
  for (unsigned n=0;n<statList.size();n++)
  {
    unsigned numNodes = statList[n].numNodes;
    lout << "level=" << n << ":  numNodes="
         << numNodes << ",totChildren="
         << statList[n].totNumChildren ;
    //              for(unsigned g=0;g<statList[n].nodeList.size();g++)
    //              {
    //                  lout << "(" << statList[n].nodeList[g].first
    //                       << "," << statList[n].nodeList[g].second
    //                       << "),";

    //              }

    lout << std::endl
         << "     avgChildren="
         << statList[n].totNumChildren/numNodes
         << ",avgchildoverlap="
         << statList[n].childoverlap/numNodes
         << ",avg margin=" << statList[n].margin/numNodes
         << ",avg coverage="
         << statList[n].coverage/numNodes
         << std::endl;
  }
  lout << std::endl;
}
bool Rtree::checkBoundingBoxes()
{
  return checkBoundingBoxesInternal(_data->root,_data->root->box,
      0,_data,this->getSource());
}
double Rtree::calculateMargin()
{
  return computeMargin(_data->root,this->getSource());
}
double Rtree::calculateOverlap()
{
  return computeOverlap(_data->root,this->getSource());
}

bool Rtree::checkLeafLevels()
{
  return checkLeafLevelsInternal(_data->root);
}

#endif


//  Internal stuff below, used for bulk loading
class SortByCenter
{
public:
  SortByCenter(int dim) : _dim(dim) {}
  int operator()(RtreeLeafEntry const *p,
                 RtreeLeafEntry const *q) const
  {
    float p_center = p->box.min(_dim) + p->box.max(_dim);
    float q_center = q->box.min(_dim) + q->box.max(_dim);
    return p_center < q_center;
  }
private:
  int _dim;
};



// TGS stuff

// Development notes:
//    The code immediately following is a modified version of the
//    'Topdown Greedy Split' bulk loader described by Garcia et. al
//    in 'A Greedy Algorithm for Bulk Loading R-trees'.  I originally
//    implemented their algorithm as written but found that the load
//    time was prohibitive.  I thus made the following modifications:
//
//    * Instead of checking all three dimensions for optimal split,
//    use the dimension of the longest side of the MBR of all current
//    entries.
//    * Instead of checking all possible split points to find the
//    optimal one in terms of some cost function, always use the
//    middle-most split
//    * In addition, instead of completely recomputing the MBR of the
//    right and left sides after the split, approximate it using the
//    actual coordinate value at the split.
//
//    These changes reduced the load time quite significantly with
//    no measurable impact on the query times (using the surface wrapper
//    intersection checker as a benchmark).

/*! Basic split step for TGS. Splits the input data into a number of
    chunks, hopefully producing reasonbly disjoint MBRS of the chunks.  
    All but the last chunk will have size M.  The total number of
    chunks is returned. */
static int TGSBasicSplit(std::vector<RtreeLeafEntry const *> &sortvector,
                         int start, int end, int M,
                         RtreeBox const &inputbox)
{
  // compute total number of elements in this split
  int n = end - start;
  if (n <= M)
    return 1; // done splitting

  CheckInterrupts::checkCurrent();

  // compute number of chunks and splits
  int num_chunks = (int)ceil((float)n / (float)M);
  int num_splits = num_chunks - 1;

  // get longest dimension of input box to use as split axis
  int dim = inputbox.getLongestDimension();

  // always use the middle-most split
  int splitnum = num_splits / 2 + 1;

  // find actual split point in vector
  int split = start + splitnum*M;

  // sort based on center values in split axis direction
  // note, actual sorting is overkill here, using nth_element instead
  // is sufficient and much faster
  std::nth_element(sortvector.begin()+start,
                   sortvector.begin()+split,
                   sortvector.begin()+end,
                   SortByCenter(dim));

  // find coordinate value at split location
  float split_location = 0.5*(sortvector[split]->box.min(dim) +
      sortvector[split]->box.max(dim));

  // approximate left and right MBRs for recursion.  Actual MBRs may
  // be smaller, but can't be bigger.
  RtreeBox left(inputbox);
  left.max(dim) = split_location;

  RtreeBox right(inputbox);
  right.min(dim) = split_location;

  // recurse to left and right sides of split until desired size is reached
  TGSBasicSplit(sortvector, start, split, M, left);
  TGSBasicSplit(sortvector, split, end, M, right);

  return num_chunks;
}

/*! Wrapper around above basic split -- computes MBR of input entries
    before invoking split. */
static int TGSBasicSplit(std::vector<RtreeLeafEntry const *> &sortvector,
                         int start, int end, int M)
{
  RtreeBox box;
  for (int k=start; k < end; ++k)
    box.expand(sortvector[k]->box);
  return TGSBasicSplit(sortvector, start, end, M, box);
}


/*! Main split step.  Splits the data up using TGSBasicSplit and
    builds appropriate R-tree nodal structure. */
void TGSSplit(std::vector<RtreeLeafEntry const *> &sortvector,
              int start,
              int end,
              int M,
              int N,
              RtreeNonLeafNode *node)
{
  CheckInterrupts::checkCurrent();

  // call basic split and compute number of chunks
  int num_chunks = TGSBasicSplit(sortvector, start, end, M);

  // compute number of entries in a fully packed child node
  int child_M = M / N;

  if (child_M == 1)
  {
    // children are leaves -- populate children and return
    for (int i=0; i<num_chunks; ++i)
    {
      int child_start = start + M*i;
      int child_end = child_start + M;
      if (child_end > end)
        child_end = end;
      RtreeLeafNode *child = new RtreeLeafNode;
      for (int j=child_start; j < child_end; ++j)
        child->addChild(*sortvector[j]);
      node->addChild(child);
    }
  }
  else
  {
    // children are non-leaves -- recurse
    for (int i=0; i<num_chunks; ++i)
    {
      RtreeNonLeafNode *child = new RtreeNonLeafNode;
      int child_start = start + M*i;
      int child_end = child_start + M;
      if (child_end > end)
        child_end = end;
      TGSSplit(sortvector, child_start, child_end, child_M, N,
               child);
      node->addChild(child);
    }
  }
}

/*! Invoke the TGS bulk loader based in given entries and node size. */
static RtreeNode *TGSBulkLoad(std::vector<RtreeLeafEntry> const &entries,
                              int N)
{
  std::vector<RtreeLeafEntry const *> sortvector(entries.size());
  for (unsigned i = 0; i<entries.size(); ++i)
    sortvector[i] = &entries[i];

  // Specific conversion required since AIX C++ couldn't find best match
  // for log(N) and pow(n, leaf_depth).  WRO 2004-Dec-21
  int leaf_depth = (int)ceil(log((double)entries.size())/log((double)N)) - 1;

  RtreeNonLeafNode *root = new RtreeNonLeafNode;

  // M is the number of entries in a packed child at current depth
  int M = (int)pow((double)N, leaf_depth);

  TGSSplit(sortvector, 0, entries.size(), M, N, root);

  return root;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
/*
    The functions below are used to implement insertion and deletion of 
    single items from an Rtree.

    The insertion algorithm is based upon that listed in the paper
    in "The R*-tree:  An Efiicient and Robust Access Method for Points and 
    Rectangles" by N. Beckmann, H. Kriegel, R. Schneider and B. Seeger
    ACM 1990

    The deletion algorithm is based upon an earlier paper:
    "R-trees:  A Dynamic Index Structure for Spatial Searching"
    by Antonin Guttman  ACM 1984

 */
static RtreeNode *insertAtLevel(RtreeNode *root,RtreeNode* newNode,
                                unsigned level,
                                unsigned *reInsertLevel,
                                Rtree::Data *data,
                                SpatialIndexSource const &source);

static RtreeNode *insertEntry(RtreeNode *root,RtreeLeafEntry &newEntry,
                              unsigned *reinsertLevel,
                              Rtree::Data *data,
                              SpatialIndexSource const &source);

typedef std::vector<RtreeNode*> RtreeNodePath;
typedef std::vector<unsigned> SortIndices,AxisSort;
typedef std::vector<AxisSort> AxesSortList;

//used for sorting a particular axis according to is lower coordinate
class AxisSortLower {
public:
  AxisSortLower(int axis,RtreeNode *parent,
                SpatialIndexSource const &source):
                  _axis(axis),_parent(parent),_source(source)
  {};
  int operator()(unsigned index1,unsigned index2)
  {
    RtreeBox box1 = _parent->getChildBox(index1);
    RtreeBox box2 = _parent->getChildBox(index2);
    return box1.getMinAxisValue(_axis) <
        box2.getMinAxisValue(_axis);
  }

private:
  int _axis;
  RtreeNode *_parent;
  SpatialIndexSource const &_source;
};
//used for sorting a particular axis according to is upper coordinate
class AxisSortUpper {
public:
  AxisSortUpper(int axis,RtreeNode *parent,
                SpatialIndexSource const &source):
                  _axis(axis),_parent(parent),_source(source){};
  int operator()(unsigned index1,unsigned index2)
  {
    RtreeBox box1 = _parent->getChildBox(index1);
    RtreeBox box2 = _parent->getChildBox(index2);
    return box1.getMaxAxisValue(_axis) <
        box2.getMaxAxisValue(_axis);

  }

private:
  int _axis;
  RtreeNode *_parent;
  SpatialIndexSource const &_source;
};

/*
    class used for sorting boxes based upon center distances
    which are placed in _centerDists for each child node
 */
class CenterSort {
public:
  CenterSort(std::vector<double> const &centerDists):
    _centerDists(centerDists)
  {};
  int operator()(unsigned index1,unsigned index2)
  {
    return _centerDists[index1] >
    _centerDists[index2];
  }

private:
  std::vector<double> const &_centerDists;
};

/*
    sorts the children of this node twice along each axis -
    first using the lower coordinate value for that axis
    second using the upper coordinate value for that axis
    in 3d this creates 6 sorts.  The sorting is just sorting
    the indices into the children rather than the children 
    themselves
 */
static void sortAlongAxes(RtreeNode *parent,
                          AxesSortList &lowerSorts,
                          AxesSortList &upperSorts,
                          SpatialIndexSource const &source)
{
  AxisSort emptySort;
  unsigned numDims = parent->getBox().getNumDimensions();
  for (unsigned axis=0;axis < numDims;axis++)
  {
    lowerSorts.push_back(emptySort);
    upperSorts.push_back(emptySort);
    for (unsigned k=0;k<parent->getNumChildren();k++)
    {
      lowerSorts[axis].push_back(k);
      upperSorts[axis].push_back(k);
    }
    std::sort(lowerSorts[axis].begin(),lowerSorts[axis].end(),
              AxisSortLower(axis,parent,source));
    std::sort(upperSorts[axis].begin(),upperSorts[axis].end(),
              AxisSortUpper(axis,parent,source));
  }
}
/*
    childIndices specifies a sorted order of this distribution
    childIndices[0] gives an index of the first child node in
    the distribution.
    So we take parent->getChildBox(childIndices[0]) to 
    parent->getChildBox(childIndices[dist-1]) as the first distribution
    and parent->getChildBox(childIndices[dist]) to 
    parent->getChildBox(childIndices[size-1]) as the other distribution
    We are computing the bouding boxes of each of these distributions
    and will use these bboxes to determine if this is the distribution we
    want to use when we split the node
 */
static void computeDistBoxes(RtreeNode *parent,unsigned dist,
                             AxisSort &childIndices,
                             RtreeBox &box1, RtreeBox &box2,
                             SpatialIndexSource const &source)
{
  box1.clear();
  box2.clear();
  for (unsigned i=0;i<dist;i++)
  {
    box1.expand(parent->getChildBox(childIndices[i]));
  }
  for (unsigned i=dist;i < childIndices.size();i++)
  {
    box2.expand(parent->getChildBox(childIndices[i]));
  }
}
// the margin is just the sum of the lengths of the sides of the bbox
static double computeDistMargin(RtreeNode *parent,unsigned dist,
                                AxisSort &childIndices,
                                SpatialIndexSource const &source)
{
  RtreeBox box1,box2;
  computeDistBoxes(parent,dist,childIndices,box1,box2,source);
  return box1.computeMargin() + box2.computeMargin();
}

static unsigned chooseSplitAxis(RtreeNode *parent,
                                AxesSortList &lowerSorts,
                                AxesSortList &upperSorts,
                                SpatialIndexSource const &source,
                                Rtree::Data *data)
{
  double minS = std::numeric_limits<double>::max();
  unsigned minAxis=0;
  unsigned numDims = parent->getBox().getNumDimensions();
  for (unsigned axis=0;axis < numDims;axis++)
  {
    double S=0.0;
    for (unsigned k=data->minN;k<=data->maxN+1-data->minN;k++)
    {
      S+=computeDistMargin(parent,k,lowerSorts[axis],source);
      S+=computeDistMargin(parent,k,upperSorts[axis],source);
    }
    if (S < minS)
    {
      minS = S;
      minAxis = axis;
    }
  }
  return minAxis;

}
/*
    computes the coverage of the given distribution.  The distribution is
    1 to dist and dist to size using childIndices as the ordering
 */
static double computeDistCoverage(RtreeNode *parent,unsigned dist,
                                  AxisSort &childIndices,
                                  SpatialIndexSource const &source)
{
  RtreeBox box1,box2;
  computeDistBoxes(parent,dist,childIndices,box1,box2,source);
  return box1.coverage() + box2.coverage();
}
/* computes the overlap between two distributions of nodes
     i.e. from 1 to dist and from dist to maxN+1 - it uses childIndices
     as the sorted indices into the children to determine the dist
 */
static double computeDistOverlap(RtreeNode *parent,unsigned dist,
                                 AxisSort &childIndices,
                                 SpatialIndexSource const &source)
{
  RtreeBox box1,box2,intersectBox;
  computeDistBoxes(parent,dist,childIndices,box1,box2,source);
  if (box1.computeIntersection(box2,intersectBox))
    return intersectBox.coverage();
  else
    return 0.0;
}
/*
    chooses a split based upon the smallest overlap
 */
static void chooseSplitIndex(RtreeNode *parent,
                             AxisSort &lowerSort,
                             AxisSort &upperSort,
                             unsigned &bestIndex,
                             AxisSort **bestSort,
                             SpatialIndexSource const &source,
                             Rtree::Data *data)
{
  AxisSort *sorts[2]={&lowerSort,&upperSort};
  double overlap,coverage,minOverlap,minCoverage=-1.0;
  minOverlap = std::numeric_limits<double>::max();
  for (unsigned k=data->minN;k<=data->maxN+1-data->minN;k++)
  {
    for (unsigned s=0;s<2;s++)
    {
      overlap=computeDistOverlap(parent,k,*(sorts[s]),source);
      if (overlap < minOverlap)
      {
        minOverlap = overlap;
        minCoverage = -1.0;
        *bestSort = sorts[s];
        bestIndex = k;
      }
      else if (overlap==minOverlap)
      {
        coverage =
            computeDistCoverage(parent,k,*(sorts[s]),source);
        if (minCoverage < 0.0)
        {
          minCoverage =
              computeDistCoverage(parent,bestIndex,
                  **bestSort,source);
        }
        if (coverage < minCoverage)
        {
          minCoverage = coverage;
          *bestSort = sorts[s];
          bestIndex = k;
        }
      }
    }
  }
}
/*
    splits a leaf node similar to splitNonLeafNode - the split criteria are
    slightly different for leaves.
 */
static RtreeNode *splitLeafNode(RtreeNode *leaf,
                                Rtree::Data *data,
                                SpatialIndexSource const &source)
{
  unsigned axis = 0, bestIndex = 0;
  AxesSortList lowerSorts,upperSorts;

  AxisSort *bestSort = 0;
  sortAlongAxes(leaf,lowerSorts,upperSorts,source);
  axis = chooseSplitAxis(leaf,lowerSorts,upperSorts,source,data);
  chooseSplitIndex(leaf,lowerSorts[axis],upperSorts[axis],
                   bestIndex,&bestSort,source,data);

  RtreeLeafNode *leafWithBox = assert_cast<RtreeLeafNode*>(leaf);
  RtreeLeafNode *newNode = new RtreeLeafNode(bestIndex+1);
  std::vector<RtreeLeafEntry> tempChildren;
  for (unsigned i=0;i<bestIndex;i++)
  {
    RtreeLeafEntry const &childEntry =
        leafWithBox->getChild((*bestSort)[i]);
    newNode->addChild(childEntry);
  }
  for (unsigned i=bestIndex;i<bestSort->size();i++)
  {
    tempChildren.push_back(leafWithBox->getChild((*bestSort)[i]));
  }
  leafWithBox->clearChildren();
  if (tempChildren.size()==0)
    throw Exception("RtreeInsDelAlg::can't"
        " have empty leaf node");
  for (unsigned i=0;i<tempChildren.size();i++)
  {
    leafWithBox->addChild(tempChildren[i]);
  }
  return newNode;
}

/*
    splits the nonleaf node by first 
    sorting along each axis the children according to their upper coordinate 
    and lower coordinate so we end up with two sorts for each axis or 6 sorts
    in 3d. 
    Then chooseSplitAxis chooses the axis with the smallest Margin (i.e. the
    sort along which we have the squarest boxes

    then it chooses which to split by finding the split that creates the
    least overlap

 */
static RtreeNonLeafNode *splitNonLeafNode(RtreeNonLeafNode *parent,
                                          Rtree::Data *data,
                                          SpatialIndexSource const &source)
{
  unsigned axis = 0, bestIndex = 0;
  AxesSortList lowerSorts,upperSorts;
  std::vector<RtreeNode*> tempChildren;
  AxisSort *bestSort = 0;
  sortAlongAxes(parent,lowerSorts,upperSorts,source);
  axis = chooseSplitAxis(parent,lowerSorts,upperSorts,source,data);
  chooseSplitIndex(parent,lowerSorts[axis],upperSorts[axis],
                   bestIndex,&bestSort,source,data);
  RtreeNonLeafNode *newNode = new RtreeNonLeafNode(bestIndex+1);

  for (unsigned i=0;i<bestIndex;i++)
  {
    newNode->addChild(parent->getChild((*bestSort)[i]));
  }
  for (unsigned i=bestIndex;i<bestSort->size();i++)
  {
    tempChildren.push_back(parent->getChild((*bestSort)[i]));
  }
  parent->clearChildren();
  if (tempChildren.size()==0)
    throw Exception("RtreeInsDelAlg::can't"
        " have node in Rtree");

  for (unsigned i=0;i<tempChildren.size();i++)
  {
    parent->addChild(tempChildren[i]);
  }

  return newNode;
}
/*
    called by overflowtreatment to split the node when it has
    more children than the maximum allowed number of children
 */
static RtreeNode *splitNode(RtreeNode *parent,
                            Rtree::Data *data,
                            SpatialIndexSource const &source)
{
  if (parent->isLeaf())
    return splitLeafNode(parent,data,source);
  else
    return splitNonLeafNode(assert_cast<RtreeNonLeafNode*>(parent),data,
        source);
}
/*
    sorts the indices of the children of node so the child whose center
    is farthest from the center of the parent is first
 */
static void sortForReinsert(RtreeNode *node,
                            SortIndices &sortIndices,
                            SpatialIndexSource const &source)
{
  std::vector<double> centerDistances;

  RtreeBox box = node->getBox();
  unsigned int const maxDims = 3;
  unsigned int const numDims = box.getNumDimensions();
  float center[maxDims];
  float center2[maxDims];
  box.getCenter(center);
  sortIndices.clear();
  for (unsigned i=0;i<node->getNumChildren();i++)
  {
    double dist,distSq;
    distSq = 0.0;
    RtreeBox const &childBox = node->getChildBox(i);
    childBox.getCenter(center2);
    for (unsigned j=0;j<numDims;j++)
    {
      dist = center[j]-center2[j];
      distSq += dist*dist;
    }
    centerDistances.push_back(distSq);
    sortIndices.push_back(i);
  }

  std::sort(sortIndices.begin(),sortIndices.end(),
            CenterSort(centerDistances));


}

static RtreeNode* reinsertLeafChildren(RtreeNode *root,
                                       RtreeNode *node,
                                       unsigned *reinsertLevel,
                                       RtreeNodePath &subtreePath,
                                       SortIndices const &sortIndices,
                                       Rtree::Data *data,
                                       SpatialIndexSource const &source)
{
  RtreeLeafEntry entry;
  std::vector<RtreeLeafEntry> reinsertList;

  for (unsigned i=0;i<node->getNumChildren();i++)
  {
    entry = assert_cast<RtreeLeafNode*>(node)->getChild(sortIndices[i]);
    reinsertList.push_back(entry);
  }
  node->clearChildren();
  for (unsigned n=data->reInsP;n<reinsertList.size();n++)
  {
    assert_cast<RtreeLeafNode*>(node)->addChild(reinsertList[n]);
  }
  for (unsigned i=subtreePath.size(); i>0;i--)
  {
    node = subtreePath[i-1];
    assert_cast<RtreeNonLeafNode*>(node)->updateBox();
  }
  for (unsigned i=data->reInsP;i>0;i--)
  {
    root = insertEntry(root,reinsertList[i-1],
        reinsertLevel,data,source);
  }
  return root;
}

/*
    removes a certain number of entries from this node and
    reinserts them at the same level into the tree

 */
static RtreeNode* reinsertNonLeafChildren(RtreeNode *root,
                                          RtreeNonLeafNode *node,
                                          unsigned nodeLevel,
                                          unsigned *reinsertLevel,
                                          RtreeNodePath &subtreePath,
                                          SortIndices const &sortIndices,
                                          Rtree::Data *data,
                                          SpatialIndexSource const &source)
{
  RtreeNode* entry;
  std::vector<RtreeNode*> reinsertList;

  for (unsigned i=0;i<node->getNumChildren();i++)
  {
    entry = node->getChild(sortIndices[i]);
    reinsertList.push_back(entry);
  }
  node->clearChildren();
  for (unsigned n=data->reInsP;n<reinsertList.size();n++)
  {
    node->addChild(reinsertList[n]);
  }
  for (unsigned i=subtreePath.size(); i>0;i--)
  {
    node = assert_cast<RtreeNonLeafNode*>(subtreePath[i-1]);
    node->updateBox();
  }
  for (unsigned i=data->reInsP;i>0;i--)
  {
    RtreeNode *prevRoot = root;
    root = insertAtLevel(root,reinsertList[i-1],nodeLevel+1,
        reinsertLevel,data,source);
    if (prevRoot != root)
      nodeLevel+=1; //just grew the tree by 1
  }
  return root;
}
/*
    removes a certain number of children from this tree and
    reinserts them at the same level

    first it sorts the child nodes based upon the distance of their centers from
    the center of the parent node.  Then it reinserts the data->reinsP 
    children furthest from the center - i.e. an attempt to shrink this nodes
    coverage - it reinserts beginning with the closest far node
 */
static RtreeNode *reinsertChildren(RtreeNode *root,
                                   RtreeNode *node,
                                   unsigned nodeLevel,
                                   unsigned *reinsertLevel,
                                   RtreeNodePath &subtreePath,
                                   Rtree::Data *data,
                                   SpatialIndexSource const &source)
{
  SortIndices sortedIndices;
  sortForReinsert(node,sortedIndices,source);

  if (node->isLeaf())
    root = reinsertLeafChildren(root,node,reinsertLevel,subtreePath,
        sortedIndices,data,source);
  else
    root = reinsertNonLeafChildren(root,
        assert_cast<RtreeNonLeafNode*>(node),
        nodeLevel,
        reinsertLevel,subtreePath,
        sortedIndices,data,source);
  return root;
}
/*
    called if insertion into node causes the number of children to exceed the
    maximum.  If so then some children are reinserted unless this function
    has already been called at this level in the tree, then it splits the node
    see paper referenced above
 */
static void overflowTreatment(RtreeNode **root,
                              RtreeNode **newNode,
                              RtreeNode *node,
                              unsigned nodeLevel,
                              unsigned *reinsertLevel,
                              RtreeNodePath &subtreePath,
                              Rtree::Data *data,
                              SpatialIndexSource const &source)
{
  *newNode = (RtreeNode*)0;
  if (node->getNumChildren() <= data->maxN)
    return;
  if (reinsertLevel)
  {
    if (nodeLevel && (*reinsertLevel > nodeLevel))
    {
      *reinsertLevel = *reinsertLevel - 1;
      *root = reinsertChildren(*root,node,nodeLevel,reinsertLevel,
          subtreePath,data,source);
    }
    else
      *newNode = splitNode(node,data,source);
  }
  else
  {
    unsigned level = nodeLevel;
    *root = reinsertChildren(*root,node,nodeLevel,&level,
        subtreePath,data,source);
  }
}
/*
    called to adjust the tree after inserting a node into the tree
    called recursively when the insertion will cause a node to be split
 */
static RtreeNode *adjustTree(RtreeNode *root,
                             unsigned *reinsertLevel,
                             RtreeNodePath &subtreePath,
                             Rtree::Data *data,
                             SpatialIndexSource const &source)
{

  if (subtreePath.empty())
  {
    return root;
  }
  RtreeNode *newNode=0;
  RtreeNode *checkNode = subtreePath.back();
  subtreePath.pop_back();
  if (checkNode->getNumChildren() > data->maxN)
  {
    overflowTreatment(&root,&newNode,checkNode,subtreePath.size(),
        reinsertLevel,subtreePath,data,source);

  }
  if (newNode)
  {
    if (root==checkNode)
    {
      RtreeNonLeafNode* newRoot = new RtreeNonLeafNode(2);
      newRoot->addChild(checkNode);
      newRoot->addChild(newNode);
      if (reinsertLevel)
        (*reinsertLevel) += 1;//just grew the tree so levels change
      return newRoot;
    }
    else
    {
      RtreeNode *parent = subtreePath.back();
      assert_cast<RtreeNonLeafNode*>(parent)->addChild(newNode);
      assert_cast<RtreeNonLeafNode*>(parent)->updateBox();
      root = adjustTree(root,reinsertLevel,
          subtreePath,data,source);

    }
  }
  else
  {
    //could recurse but why, the rest of the parents just
    //need box updated
    //if we reinserted elsewhere that is ok we might
    //be updating parents that don't need updating
    //but sometimes they will need updating
    //all of these will be non leaf nodes as we popped the back
    //before getting here
    while (subtreePath.size())
    {
      RtreeNode* parent=subtreePath.back();
      assert_cast<RtreeNonLeafNode*>(parent)->updateBox();
      subtreePath.pop_back();
    }
  }
  return root;
}

/*
    returns the extra spatial coverage (volume in 3d) the node will use when
    it is expanded to include newBox
 */
static double computeExpansionCost(RtreeNode* child,RtreeBox const &newBox)
{
  double oldCoverage,newCoverage;
  double expansion;
  RtreeBox tempBox = child->getBox();
  oldCoverage = tempBox.coverage();
  tempBox.expand(newBox);
  newCoverage = tempBox.coverage();
  expansion = newCoverage-oldCoverage;
  return expansion;
}
/*
    computes how much the child with the given index will overlap with
    its siblings if it is expanded to include newBox

    This overlap is computed as the space occupied (volume for 3d boxes)
    by intersection of the box with each of its siblings
 */
static double computeLeafExpansionOverlap(RtreeNonLeafNode *node,
                                          unsigned leafIndex,
                                          RtreeBox const &newBox)
{
  RtreeBox overlapBox;
  RtreeNode *child2,
  *leaf = node->getChild(leafIndex);
  RtreeBox box1 = leaf->getBox();
  double cost = 0.0;
  unsigned numChildren = node->getNumChildren();
  box1.expand(newBox);
  for (unsigned j=0;j<numChildren;j++)
  {
    if (j==leafIndex)
      continue;
    child2 = node->getChild(j);
    RtreeBox box2 = child2->getBox();
    if (box2.computeIntersection(box1,overlapBox))
      cost+=overlapBox.coverage();
  }
  return cost;
}

/*
    using newBox finds the best child of node in which newBox would be
    inserted.  This is used to choose the insertion path for inserting
    an entry with a bbox newBox into the tree

    It bases this for Leaf Nodes on how much overlap adding the box
    to a particular child will cause with its other siblings.  The
    sibling with the lowest overlap is used

    for non-leaf nodes it chooses the child that will need the smallest
    amount of expansion to include the new box

    returns the preferred child

    See paper listed above
 */
static RtreeNode *findBestChildForInsert(RtreeNonLeafNode *node,
                                         RtreeBox &newBox)
{
  double minCost,cost,minExpand=-1.0,minCoverage=-1.0,expand,coverage;
  RtreeNode *bestChild=0,*child1;
  unsigned numChildren = node->getNumChildren();
  minExpand=minCost = std::numeric_limits<double>::max();
  child1 = node->getChild(0);
  if (child1->isLeaf())
  {
    //for leaves find child whose new box makes
    //the smallest amount of overlap with siblings
    for (unsigned i=0;i<numChildren;i++)
    {
      //this is O(numChildren^2) paper discusses possibility
      //of approximating this if numChildren becomes large
      //but I am not going to do the approximation right now
      child1 = node->getChild(i);
      cost = computeLeafExpansionOverlap(node,i,newBox);
      if (cost < minCost)
      {
        minCost = cost;
        bestChild = child1;
        minExpand = -1.0;
      }
      else if (cost == minCost)
      {
        if (minExpand < 0.0)
          minExpand = computeExpansionCost(bestChild,newBox);
        expand = computeExpansionCost(child1,newBox);
        if (expand < minExpand)
        {
          minExpand = expand;
          bestChild = child1;
          minCoverage = -1.0;
        }
        else if (expand == minExpand)
        {
          if (minCoverage < 0.0)
          {
            RtreeBox childBox = bestChild->getBox();
            childBox.expand(newBox);
            minCoverage = childBox.coverage();
          }
          RtreeBox childBox = child1->getBox();
          childBox.expand(newBox);
          coverage = childBox.coverage();
          if (coverage < minCoverage)
          {
            minCoverage = coverage;
            bestChild = child1;
          }
        }
      }
    }
  }
  else
  {
    //for non-leaves find child whose expansion cost is smallest
    //i.e. the additional space covered by adding child is smallest
    for (unsigned i=0;i<numChildren;i++)
    {
      child1 = node->getChild(i);
      cost = computeExpansionCost(child1,newBox);
      if (cost < minCost)
      {
        minCoverage = -1.0;
        minCost = cost;
        bestChild = child1;
      }
      else if (cost==minCost)
      {
        if (minCoverage < 0.0)
        {
          RtreeBox childBox = bestChild->getBox();
          childBox.expand(newBox);
          minCoverage = childBox.coverage();
        }
        RtreeBox childBox = child1->getBox();
        childBox.expand(newBox);
        coverage = childBox.coverage();
        if (coverage < minCoverage)
        {
          minCoverage = coverage;
          bestChild = child1;
        }
      }
    }
  }
  return bestChild;
}
/*
    chooses the best leaf in which to insert newBox
    returns the path from the root down to newBox
 */
static RtreeNode *chooseLeaf(RtreeNode *node, RtreeBox &newBox,
                             RtreeNodePath &subtreePath)
{
  RtreeNode *bestChild;

  subtreePath.push_back(node);
  if (node->isLeaf())
    return node;

  bestChild = findBestChildForInsert(assert_cast<RtreeNonLeafNode*>(node),
      newBox);

  RtreeNode *newNode = chooseLeaf(bestChild,newBox,subtreePath);
  return newNode;
}

/*
    searches for the entry in delEntry.  returns true if found false otherwise
    upon returning true, entryOwner will be the parent node that owns the entry
    entryIndex the index of entry into the entryOwner and subtreePath will
    contain the path from top of the tree to the leaf
 */
static bool findEntry(RtreeNode *node, RtreeLeafEntry &delEntry,
                      RtreeNode **entryOwner,
                      unsigned &entryIndex,
                      RtreeNodePath &subtreePath,
                      SpatialIndexSource const &source,
                      bool (RtreeBox::*compare)(RtreeBox const &) const)
{
  subtreePath.push_back(node);
  if (node->isLeaf())
  {
    RtreeLeafNode *leafNode = assert_cast<RtreeLeafNode*>(node);
    if (leafNode->findEntry(delEntry.id,entryIndex))
    {
      *entryOwner = leafNode;
      return true;
    }
  }
  else
  {
    for (unsigned i=0;i<node->getNumChildren();i++)
    {
      RtreeBox const &box = node->getChildBox(i);
      if ((box.*compare)(delEntry.box))
      {
        if (findEntry(assert_cast<RtreeNonLeafNode*>(node)->getChild(i),
            delEntry,
            entryOwner,entryIndex,subtreePath,source,compare))
          return true;
      }
    }
  }

  subtreePath.pop_back();
  return false;
}

/*
    searches for the entry in delEntry.  returns true if found false otherwise
    upon returning true, entryOwner will be the parent node that owns the entry
    entryIndex the index of entry into the entryOwner and subtreePath will
    contain the path from top of the tree to the leaf
 */
static bool findInvalidEntry(RtreeNode *node, RtreeLeafEntry &delEntry,
                             RtreeNode **entryOwner,
                             unsigned &entryIndex,
                             RtreeNodePath &subtreePath)
{
  subtreePath.push_back(node);
  if (node->isLeaf())
  {
    RtreeLeafNode *leafNode = assert_cast<RtreeLeafNode*>(node);
    if (leafNode->findEntry(delEntry.id,entryIndex))
    {
      *entryOwner = leafNode;
      return true;
    }
  }
  else
  {
    for (unsigned i=0;i<node->getNumChildren();i++)
    {
      if (findInvalidEntry(assert_cast<RtreeNonLeafNode*>(node)->getChild(i),
          delEntry,
          entryOwner,
          entryIndex,
          subtreePath))
        return true;
    }
  }

  subtreePath.pop_back();
  return false;
}


static bool findEntry(RtreeNode *node, RtreeLeafEntry &delEntry,
                      RtreeNode **entryOwner,
                      unsigned &entryIndex,
                      RtreeNodePath &subtreePath,
                      SpatialIndexSource const &source)
{
  // first search using strict containment
  if (findEntry(node, delEntry, entryOwner, entryIndex,
      subtreePath, source, &RtreeBox::contains))
    return true;
  // not found so search all nodes with intersecting boxes
  if (findEntry(node, delEntry, entryOwner, entryIndex,
      subtreePath, source, &RtreeBox::intersects))
    return true;
  // still not found - call brute force method
  if (findInvalidEntry(node, delEntry, entryOwner, entryIndex,
      subtreePath))
    return true;
  return false;
}


/*
    chooses the best subtree for inserting the box at the given level
    similar to chooseLeaf but allows insertion at higher levels in the tree
 */
static RtreeNode *chooseSubtree(RtreeNode *node, RtreeBox &newBox,
                                RtreeNodePath &subtreePath,unsigned level)
{
  RtreeNode *bestChild;

  subtreePath.push_back(node);

  if ((level < subtreePath.size()) ||
      (node->isLeaf()))
    return node;

  bestChild = findBestChildForInsert(assert_cast<RtreeNonLeafNode*>(node),
      newBox);

  RtreeNode *newNode =
      chooseSubtree(bestChild,newBox,subtreePath,level);
  return newNode;
}
/*
    inserts the node at the given level.  No checks are made to see if the
    insertion will maintain a balanced tree it is assumed that the caller
    knows that it will
 */
static RtreeNode *insertAtLevel(RtreeNode *root,RtreeNode* newNode,
                                unsigned level,
                                unsigned *reInsertLevel,
                                Rtree::Data *data,
                                SpatialIndexSource const &source)
{
  RtreeNodePath subtreePath;
  RtreeBox newBox;
  RtreeNode *node = chooseSubtree(root,newNode->box,subtreePath,level-1);
  if (node->isLeaf())
    throw Exception("RtreeInsDelAlg:: should never call"
        "  on leaf node");
  assert_cast<RtreeNonLeafNode*>(node)->addChild(newNode);

  root = adjustTree(root,reInsertLevel,
      subtreePath,data,source);

  return root;
}

static RtreeNode *condenseTree(RtreeNodePath &subtreePath,
                               Rtree::Data *data,
                               SpatialIndexSource const &source)
{
  unsigned index;
  std::vector<RtreeNode*> orphanedNodes;
  RtreeNode *parent = subtreePath.back();
  RtreeNode* root = subtreePath.front();
  subtreePath.pop_back();

  while ((parent != root) &&
      (parent->getNumChildren() < data->minN))
  {
    RtreeNonLeafNode *newParent =
        assert_cast<RtreeNonLeafNode*>(subtreePath.back());
    if (newParent->findChild(parent,index))
      newParent->removeChild(index);
    else
      throw Exception("RtreeInsDelAlg::condenseTree should "
          "always find child node");
    orphanedNodes.push_back(parent);
    parent = newParent;
    subtreePath.pop_back();

  }

  unsigned insertLevel = subtreePath.size() + 1;

  while (subtreePath.size())
  {
    RtreeNonLeafNode* nonLeafParent =
        assert_cast<RtreeNonLeafNode*>(subtreePath.back());
    nonLeafParent->updateBox();
    subtreePath.pop_back();
  }

  if (orphanedNodes.size())
  {
    for (unsigned i=orphanedNodes.size()-1;i>0;i--)
    {
      insertLevel+=1;
      for (unsigned j=0;j<orphanedNodes[i]->getNumChildren();j++)
      {
        RtreeNode* prevRoot = root;
        RtreeNode* nodeToInsert =
            assert_cast<RtreeNonLeafNode*>(orphanedNodes[i])->getChild(j);
        root = insertAtLevel(root,nodeToInsert,
            insertLevel,0,data,source);
        if (prevRoot != root)
        {
          //tree grew by one level so need to increment insertLevel
          insertLevel+=1;
        }
      }
      orphanedNodes[i]->clearChildren();

      delete orphanedNodes[i];
    }
    RtreeNode *leaf = orphanedNodes[0];
    if (!leaf->isLeaf())
      throw Exception("RtreeInsDelAlg:: should always"
          "  be a leaf node");

    RtreeLeafEntry entry;
    for (unsigned i=0;i<leaf->getNumChildren();i++)
    {
      entry = assert_cast<RtreeLeafNode*>(leaf)->getChild(i);
      root = insertEntry(root,entry,0,data,source);
    }
    leaf->clearChildren();
    delete leaf;
  }
  if (!root->isLeaf() && (root->getNumChildren()==1))
  {
    RtreeNode *newRoot = assert_cast<RtreeNonLeafNode*>(root)->getChild(0);
    root->clearChildren();
    delete root;
    root = newRoot;
  }
  return root;
}

/*
    This is only called by Rtree::deleteItem

    deletes the RtreeLeafEntry from the tree whose root is root

    root - root of the Rtree
    delEntry - entry to be deleted id and bbox
    data - pointer to Rtree::Data field for minN and maxN
    source - spatial index source
 */
static bool removeEntry(RtreeNode *&root,
                        RtreeLeafEntry &delEntry,
                        Rtree::Data *data,
                        SpatialIndexSource const &source,
                        bool assumeValid)
{
  RtreeNodePath subtreePath;
  RtreeNode *newLeaf=0;
  unsigned index;
  if (assumeValid)
  {
    if (!findEntry(root,delEntry,&newLeaf,index,subtreePath,source))
      return false;
  }
  else
  {
    if (!findInvalidEntry(root, delEntry, &newLeaf, index, subtreePath))
      return false;
  }

  newLeaf->removeChild(index);
  root = condenseTree(subtreePath,data,source);
  return true;
}

/*
     This is only called by Rtree::insertItem

     inserts the item -    newEntry which contains(id,bbox) into the rtree
     beginning at the root specified by root.  Based on algorithm described
     in "The R*-tree:  An Efiicient and Robust Access Method for Points and 
     Rectangles" by N. Beckmann, H. Kriegel, R. Schneider and B. Seeger

     root - root of the tree
     newEntry - contains the id and bounding box of the new item for insertion
     data - a pointer to the Rtree::Data field - contains minN and maxN the
     allowed minimum and maximum number of entries in a node
     source - the spatial index source used when RtreeNoBoxLeafNodes are 
     used in the Rtree for computing bboxes

 */
static bool insertEntry(RtreeNode *&root,
                        RtreeLeafEntry &newEntry,
                        Rtree::Data *data,
                        SpatialIndexSource const &source)
{
  root = insertEntry(root,newEntry,0,data,source);
  return true;
}

static RtreeNode *insertEntry(RtreeNode *root,RtreeLeafEntry &newEntry,
                              unsigned *reinsertLevel,
                              Rtree::Data *data,
                              SpatialIndexSource const &source)
{
  RtreeNodePath subtreePath;
  RtreeBox newBox;
  RtreeNode *leaf = chooseLeaf(root,newEntry.box,subtreePath);

  assert_cast<RtreeLeafNode*>(leaf)->addChild(newEntry);

  root = adjustTree(root,reinsertLevel,subtreePath,data,source);
  return root;
}


static bool removeEntities(RtreeNode *node,
                           Subset const &deletedEntities,
                           unsigned numEntitiesToDelete,
                           unsigned numDeletedEntities)
{
  bool changed = false;
  if (node->isLeaf())
  {
    RtreeLeafNode *leafNode = assert_cast<RtreeLeafNode*>(node);
    for (unsigned i = 0; i < leafNode->getNumChildren(); )
    {
      RtreeLeafEntry const &entry = leafNode->getChild(i);
      if (deletedEntities.contains(entry.id))
      {
        changed = true;

        // found entry -- delete it
        leafNode->removeChild(i);

        numDeletedEntities++;
        if (numDeletedEntities == numEntitiesToDelete)
          break; // have deleted all entities -- nothing left to do!

        // don't increment i so we can check next entry!
      }
      else
        // increment i
        ++i;
    }
  }
  else
  {
    for (unsigned i=0;i<node->getNumChildren(); )
    {
      RtreeNode *child =
          assert_cast<RtreeNonLeafNode*>(node)->getChild(i);
      bool childChanged = removeEntities(child,
          deletedEntities,
          numEntitiesToDelete,
          numDeletedEntities);
      if (childChanged)
        changed = true;

      // see if child is now empty -- if so, delete it
      if (child->getNumChildren() == 0)
      {
        node->removeChild(i);
        delete child;
      }
      else
        ++i;

      if (numDeletedEntities == numEntitiesToDelete)
        break; // have deleted all entities -- nothing left to do!

    }
  }

  if (changed)
  {
    // Something has changed at or below this node.   Tighten it's
    // bounding box.
    node->box.clear();
    for (unsigned i = 0; i<node->getNumChildren(); ++i)
      node->box.expand(node->getChildBox(i));
  }

  return changed;
}


//         finish me later???
//   static void getOrphanedEntities(RtreeNode *node,
//                                std::vector<Id> &orphanedEntities)
//   {
//     if(node->isLeaf())
//       {


//   static void removeUnderfullNodes(RtreeNode *node,
//                                 unsigned minN,
//                                 std::vector<Id> &orphanedEntities)
//   {
//     if(!node->isLeaf())
//       {
//      RtreeNonLeafNode *n = assert_cast<RtreeNonLeafNode *>(node);
//      for(unsigned i = 0; i<node.getNumChildren();)
//        {
//          // recurse first so algorithm works bottom up.
//          removeUnderfullNodes(n->getChild(i), minN, orphanedEntities);
//          if (n->getChild(i).getNumChildren() < minN)
//            {
//              // node has underflowed. Add all of it's entities to orphaned
//              // list


//       }



static void getAllEntriesInternal(RtreeNode *root,Subset &items,
                                  Subset &duplicates)
{
  if (root->isLeaf())
  {
    Id id;
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      RtreeLeafEntry entry;
      entry = assert_cast<RtreeLeafNode*>(root)->getChild(i);
      id = entry.id;
      if (items.contains(id))
        duplicates.add(id);
      else
        items.add(id);
    }
  }
  else
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      RtreeNode* child =
          assert_cast<RtreeNonLeafNode*>(root)->getChild(i);
      getAllEntriesInternal(child,items,duplicates);
    }
}



///////////////////////useful for debugging purposes but commented out
#ifdef InsertionDebugOn

static bool checkBoundingBoxesInternal(RtreeNode *root,
                                       RtreeBox &parentBox,
                                       int level,
                                       Rtree::Data *data,
                                       SpatialIndexSource const &source)
{
  RtreeBox myBox = root->getBox();
  //    if(!parentBox.contains(myBox))
  //    {
  //        lout << "level=" << level << ":  node"
  //             << " with box=" << myBox
  //             << " not contained by parent with"
  //             << " box=" << parentBox << std::endl;
  //        return false;
  //    }
  RtreeBox childrenBox;
  for (unsigned i=0;i<root->getNumChildren();i++)
  {

    RtreeBox const &childBox = root->getChildBox(i);
    childrenBox.expand(childBox);
  }
  if (!(myBox == childrenBox))
  {
    lout << "node =" << root
         << "level=" << level
         << " with childrenbox=" << childrenBox
         << " not equal parent with box="
         << myBox
         << std::endl;
    return false;
  }

  bool res = true;
  if (!root->isLeaf())
  {
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      RtreeNode* child =
          assert_cast<RtreeNonLeafNode*>(root)->getChild(i);
      if (!checkBoundingBoxesInternal(child,myBox,level+1,data,source))
      {
        res = false;
      }
    }
  }
  return res;
}

static unsigned getLeafLevel(RtreeNode *start,unsigned level)
{
  if (start->isLeaf())
    return level;
  RtreeNode* firstChild = assert_cast<RtreeNonLeafNode*>(start)->getChild(0);
  return getLeafLevel(firstChild,level+1);
}

static bool checkAllLeafLevels(RtreeNode *root, unsigned leafLevel,
                               unsigned curLevel)
{
  if (root->isLeaf())
  {
    if (curLevel != leafLevel)
    {
      lout << "Leaf=" << root << " has level=" << curLevel
           << " level should be=" << leafLevel << std::endl;
      return false;
    }
    return true;
  }
  else
  {
    bool result = true;
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      RtreeNode* childNode =
          assert_cast<RtreeNonLeafNode*>(root)->getChild(i);
      if (!checkAllLeafLevels(childNode,leafLevel,curLevel+1))
        result = false;
    }
    return result;
  }

}
static bool checkLeafLevelsInternal(RtreeNode *root)
{
  unsigned leafLevel = getLeafLevel(root,0);

  return checkAllLeafLevels(root,leafLevel,0);
}

static double computeChildOverlap(RtreeNode *node,
                                  SpatialIndexSource const &source)
{
  RtreeBox box1,box2,intBox;
  double overlap=0.0;
  for (unsigned i=0;i<node->getNumChildren();i++)
  {
    RtreeBox const &box1 = node->getChildBox(i);
    for (unsigned j=i;j<node->getNumChildren();j++)
    {
      RtreeBox const &box2 = node->getChildBox(j);
      if (box1.computeIntersection(box2,intBox))
        overlap+=intBox.coverage();
    }
  }
  return overlap;
}

static double computeOverlap(RtreeNode *root,
                             SpatialIndexSource const &source)
{
  double overlap = computeChildOverlap(root,source);
  if (!root->isLeaf())
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      overlap +=
          computeOverlap(assert_cast<RtreeNonLeafNode*>(root)->getChild(i),
              source);
    }
  return overlap;
}
static double computeMargin(RtreeNode *root,
                            SpatialIndexSource const &source)
{
  RtreeBox box1,box2,intBox;
  double margin=root->getBox().computeMargin();

  if (!root->isLeaf())
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      margin +=
          computeMargin(assert_cast<RtreeNonLeafNode*>(root)->getChild(i),
              source);
    }

  return margin;
}
static void addToLevelStats(RtreeNode *node,LevelStats &levelStats,
                            SpatialIndexSource const &source)
{
  levelStats.nodeList.push_back(
      std::pair<RtreeNode*,unsigned>(node,node->getNumChildren()));
  levelStats.numNodes+=1;
  levelStats.totNumChildren+=node->getNumChildren();
  levelStats.childoverlap+=computeChildOverlap(node,source);
  levelStats.margin += node->getBox().computeMargin();
  levelStats.coverage += node->getBox().coverage();
}

static void computeLevelStatsInt(RtreeNode *root,LevelStatList &statList,
                                 unsigned level,
                                 SpatialIndexSource const &source)
{

  if (statList.size()<=level)
    statList.push_back(LevelStats());
  addToLevelStats(root,statList[level],source);

  if (!root->isLeaf())
  {
    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      computeLevelStatsInt(assert_cast<RtreeNonLeafNode*>(root)->getChild(i),
          statList,level+1,source);
    }
  }
}

static bool findNode(RtreeNode *root, RtreeNode *node,
                     RtreeNode **nodeOwner,
                     unsigned &nodeIndex,
                     RtreeNodePath &subtreePath)
{

  if (root->isLeaf())
  {
    return false;
  }
  //    if(root->getBox().contains(node->getBox()))
  {
    subtreePath.push_back(root);

    for (unsigned i=0;i<root->getNumChildren();i++)
    {
      RtreeNode* childNode;
      RtreeNonLeafNode *nonLeaf =
          assert_cast<RtreeNonLeafNode*>(root);
      childNode=(nonLeaf)->getChild(i);
      if (childNode== node)
      {
        subtreePath.push_back(node);
        *nodeOwner = root;
        nodeIndex = i;
        return true;
      }
      else //if(childNode->getBox().contains(node->getBox()))
      {
        if (findNode(childNode,node,
            nodeOwner,nodeIndex,subtreePath))
          return true;
      }
    }

    subtreePath.pop_back();
  }
  return false;
}
static bool findNodeLevel(RtreeNode *root, RtreeNode *node,
                          unsigned &level)
{
  RtreeNode *nodeOwner;
  RtreeNodePath subtreePath;

  if (findNode(root,node,&nodeOwner,level,subtreePath))
  {
    level = subtreePath.size()-1;
    return true;
  }
  level = 0;
  return false;
}
#endif
}
