# A summary of aerial videography

## Advantage of aerial videography than humans
 - Human can only take video in a local space, while aerial vehicle can take video from more interesting points of view
 - Easily cover a large region of interests
 - Tracking fast moving object

## Main requirement of video && what is good video
 - The target object's position, size, shoting angle, salient feature
 - Cover one object from different view, attitute.
 - The trajectory and turning of camera should be smooth 

 ## Related works
 
 ### Interaction and interface

 It is very challenging to simutanously control the aerial vehicle and the camera. Actually in practice, it often require two operators, one for controlling the drone while the other control the camera. So it is necessary to desigh some effective tool for controlling and interaction. 
 
  - _XPose_ is a representive work, in which the user can interact directly with the target object. In a typical use of XPose, the use choose the object from video displaying screen. Then the drone will autonomously fly in some predefined pattern and get some sample images. After that the user can review the images and choose the best one for further modification by dragging some other object in the image. The drone will fly to the corresponding pose. The localization and object tracking is done using ORB-SLAM

  - _Toric space_ is also impressive. Instead of using pose of the drone and camera, it use (alpha, beta, gamma), some angle with respect to the object to describe to state of drone and camera. Object position, size, vantage angle can be computed efficiently. This provide a efficient method to control the drone and camera in the image space (similar to XPose)

  - Flycam is a simple gesture controlling tool to control drone and camera. 


### Video in one single sketch

When using a drone to take video for a large region, which is often needed in city view, tourist site videography. It is convient to define just some keyframe or interest points and then a trajectory is generated and executed.

Yang and Xie from Shenzhen University have to very attractive work. In "Uncut aerial video...", a sketch is given by the user in a 2.5D map. Their algorithm select some interest region and sample some pose along the sketch. The best conbination of these sample poses are chosen by solving STSP problem. Similarly, in "creating and chaining cmaera moves", the user just need to select some interest regions and start/end pose. Candidates local moves around each interest region will be generated. The best combination of candidates is dicided by solving STSP problems.

### Smoothing trajectory and camera motion
Smooth trajectory of drone and motion of camera is important. Both in Yand and Xie's work and some other work optimize the time allocation of the trajectory and camrea motion.

### Target Tracking

- Chen jing's work track a object using some arbitrary pattern.
- AIT's work use model predictive control and consider object position, size, angle and collision in their optimization.
- skydio, tracking in forest.
 
