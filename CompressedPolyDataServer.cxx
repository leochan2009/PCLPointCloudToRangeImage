/*=========================================================================

  Program:   Open IGT Link -- Example for Tracking Data Server
  Module:    $RCSfile: $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <fstream>
#include "api/svc/codec_api.h"
#include "api/svc/codec_def.h"
#include "api/svc/codec_app_def.h"
#include "utils/BufferedData.h"
#include "utils/FileInputStream.h"
#include "api/sha1.c"
#include "igtl_header.h"
#include "igtl_video.h"
#include "igtlMessageDebugFunction.h"
#include "range_image_local.cpp"
//#include <pcl/visualization/pcl_visualizer.h>

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlVideoMessage.h"
#include "igtlServerSocket.h"
#include "igtlMultiThreader.h"
#include <math.h>

#define IGTL_IMAGE_HEADER_SIZE          72

void* ThreadFunction(void* ptr);
int   SendVideoData(igtl::Socket::Pointer& socket, igtl::VideoMessage::Pointer& videoMsg);

typedef struct {
  int   nloop;
  igtl::MutexLock::Pointer glock;
  igtl::Socket::Pointer socket;
  int   interval;
  int   stop;
} ThreadData;

std::string     polyFile = "";

int main(int argc, char* argv[])
{

  //------------------------------------------------------------
  // Parse Arguments

  if (argc != 3) // check number of arguments
    {
    // If not correct, print usage
    std::cerr << "Usage: " << argv[0] << " <port>"    << std::endl;
    std::cerr << "    <port>     : Port # (18944 in Slicer default)"   << std::endl;
    std::cerr << "    <PolyFile>     : the name of the poly data with full directory "   << std::endl;
    exit(0);
    }

  int    port     = atoi(argv[1]);
  polyFile = argv[2];
  
  pcl::PointCloud<PointXYZ> pointCloud;
  
  // Generate the data
  for (float y=-0.5f; y<=0.5f; y+=0.01f) {
    for (float z=-0.5f; z<=0.5f; z+=0.01f) {
      PointXYZ point;
      point.x = 2.0f - y;
      point.y = y;
      point.z = z;
      pointCloud.points.push_back(point);
    }
  }
  pointCloud.width = (uint32_t) pointCloud.points.size();
  pointCloud.height = 1;
  float angularResolution = (float) (  1.0f * (M_PI/180.0f));  //   1.0 degree in radians
  float maxAngleWidth     = (float) (360.0f * (M_PI/180.0f));  // 360.0 degree in radians
  float maxAngleHeight    = (float) (180.0f * (M_PI/180.0f));  // 180.0 degree in radians
  Eigen::Affine3f sensorPose = (Eigen::Affine3f)Eigen::Translation3f(0.0f, 0.0f, 0.0f);
  RangeImageLocal::CoordinateFrame coordinate_frame = RangeImageLocal::CAMERA_FRAME;
  float noiseLevel=0.00;
  float minRange = 0.0f;
  int borderSize = 1;
  
  RangeImageLocal rangeImage;
  rangeImage.createFromPointCloud(pointCloud, angularResolution, maxAngleWidth, maxAngleHeight, sensorPose, coordinate_frame, noiseLevel, minRange, borderSize);
  
  std::cout << rangeImage << "\n";
  
  igtl::ServerSocket::Pointer serverSocket;
  serverSocket = igtl::ServerSocket::New();
  int r = serverSocket->CreateServer(port);

  if (r < 0)
  {
    std::cerr << "Cannot create a server socket." << std::endl;
    exit(0);
  }

  igtl::MultiThreader::Pointer threader = igtl::MultiThreader::New();
  igtl::MutexLock::Pointer glock = igtl::MutexLock::New();
  ThreadData td;

  while (1)
    {
    //------------------------------------------------------------
    // Waiting for Connection
    int threadID = -1;
    igtl::Socket::Pointer socket;
    socket = serverSocket->WaitForConnection(1000);
    
    if (socket.IsNotNull()) // if client connected
      {
      std::cerr << "A client is connected." << std::endl;

      // Create a message buffer to receive header
      igtl::MessageHeader::Pointer headerMsg;
      headerMsg = igtl::MessageHeader::New();
      //------------------------------------------------------------
      // loop
      for (;;)
        {
        // Initialize receive buffer
        headerMsg->InitPack();

        // Receive generic header from the socket
        int rs = socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize());
        if (rs == 0)
          {
          if (threadID >= 0)
            {
            td.stop = 1;
            threader->TerminateThread(threadID);
            threadID = -1;
            }
          std::cerr << "Disconnecting the client." << std::endl;
          td.socket = NULL;  // VERY IMPORTANT. Completely remove the instance.
          socket->CloseSocket();
          break;
          }
        if (rs != headerMsg->GetPackSize())
          {
          continue;
          }

        // Deserialize the header
        headerMsg->Unpack();

        // Check data type and receive data body
        if (strcmp(headerMsg->GetDeviceType(), "STT_VIDEO") == 0)
          {
          std::cerr << "Received a STT_VIDEO message." << std::endl;
        
          igtl::StartVideoDataMessage::Pointer startVideoMsg;
          startVideoMsg = igtl::StartVideoDataMessage::New();
          startVideoMsg->SetMessageHeader(headerMsg);
          startVideoMsg->AllocatePack();
        
          socket->Receive(startVideoMsg->GetPackBodyPointer(), startVideoMsg->GetPackBodySize());
          int c = startVideoMsg->Unpack(1);
          if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
            {
              td.interval = startVideoMsg->GetResolution();
              td.glock    = glock;
              td.socket   = socket;
              td.stop     = 0;
              threadID    = threader->SpawnThread((igtl::ThreadFunctionType) &ThreadFunction, &td);
            }
          }
        else if (strcmp(headerMsg->GetDeviceType(), "STP_VIDEO") == 0)
          {
          socket->Skip(headerMsg->GetBodySizeToRead(), 0);
          std::cerr << "Received a STP_VIDEO message." << std::endl;
          if (threadID >= 0)
            {
              td.stop  = 1;
              threader->TerminateThread(threadID);
              threadID = -1;
              std::cerr << "Disconnecting the client." << std::endl;
              td.socket = NULL;  // VERY IMPORTANT. Completely remove the instance.
              socket->CloseSocket();
            }
          break;
          }
        else
          {
          std::cerr << "Receiving : " << headerMsg->GetDeviceType() << std::endl;
          socket->Skip(headerMsg->GetBodySizeToRead(), 0);
          }
        }
      }
    }
    
  //------------------------------------------------------------
  // Close connection (The example code never reaches to this section ...)
  serverSocket->CloseSocket();

}


struct EncodeFileParam {
  const char* pkcFileName;
  const char* pkcHashStr;
  EUsageType eUsageType;
  int iWidth;
  int iHeight;
  float fFrameRate;
  SliceModeEnum eSliceMode;
  bool bDenoise;
  int  iLayerNum;
  bool bLossless;
  bool bEnableLtr;
  bool bCabac;
  int iTargetBitrate;
  // unsigned short iMultipleThreadIdc;
};

void EncFileParamToParamExt (EncodeFileParam* pEncFileParam, SEncParamExt* pEnxParamExt) {
  pEnxParamExt->iUsageType       = pEncFileParam->eUsageType;
  pEnxParamExt->iPicWidth        = pEncFileParam->iWidth;
  pEnxParamExt->iPicHeight       = pEncFileParam->iHeight;
  pEnxParamExt->fMaxFrameRate    = pEncFileParam->fFrameRate;
  pEnxParamExt->iSpatialLayerNum = pEncFileParam->iLayerNum;
  
  pEnxParamExt->bEnableDenoise   = pEncFileParam->bDenoise;
  pEnxParamExt->bIsLosslessLink  = pEncFileParam->bLossless;
  pEnxParamExt->bEnableLongTermReference = pEncFileParam->bEnableLtr;
  pEnxParamExt->iEntropyCodingModeFlag   = pEncFileParam->bCabac ? 1 : 0;
  pEnxParamExt->iTargetBitrate = pEncFileParam->iTargetBitrate;
  pEnxParamExt->uiMaxNalSize = 1500;
  //pEnxParamExt->uiIntraPeriod = 1;
  pEnxParamExt->iNumRefFrame = AUTO_REF_PIC_COUNT;
  if (pEncFileParam->eSliceMode != SM_SINGLE_SLICE && pEncFileParam->eSliceMode != SM_SIZELIMITED_SLICE) //SM_DYN_SLICE don't support multi-thread now
    pEnxParamExt->iMultipleThreadIdc = 4;
  
  for (int i = 0; i < pEnxParamExt->iSpatialLayerNum; i++) {
    pEnxParamExt->sSpatialLayers[i].iVideoWidth     = pEncFileParam->iWidth;
    pEnxParamExt->sSpatialLayers[i].iVideoHeight    = pEncFileParam->iHeight;
    pEnxParamExt->sSpatialLayers[i].fFrameRate      = pEncFileParam->fFrameRate;
    pEnxParamExt->sSpatialLayers[i].iSpatialBitrate = pEncFileParam->iTargetBitrate;
    //pEnxParamExt->sSpatialLayers[i].uiProfileIdc = PRO_UNKNOWN;///< value of profile IDC (PRO_UNKNOWN for auto-detection)
    //pEnxParamExt->sSpatialLayers[i].uiLevelIdc = LEVEL_UNKNOWN;///< value of profile IDC (0 for auto-detection)
    //pEnxParamExt->sSpatialLayers[i].iDLayerQp = 0;///< value of level IDC (0 for auto-detection)
    pEnxParamExt->sSpatialLayers[i].sSliceArgument.uiSliceMode = pEncFileParam->eSliceMode;
    if (pEncFileParam->eSliceMode == SM_SIZELIMITED_SLICE) {
      pEnxParamExt->sSpatialLayers[i].sSliceArgument.uiSliceSizeConstraint = 600;
      pEnxParamExt->bUseLoadBalancing = false;
    }
  }
  pEnxParamExt->iTargetBitrate *= pEnxParamExt->iSpatialLayerNum;
  pEnxParamExt->bEnableFrameSkip = true;
}

static EncodeFileParam kFileParamArray =
{
  "res/Cisco_Absolute_Power_1280x720_30fps.yuv",
  "dfd4666f9b90d5d77647454e2a06d546adac6a7c", CAMERA_VIDEO_REAL_TIME, 1280, 720, 30.0f, SM_SIZELIMITED_SLICE, false, 1, false, false, true, 5000000
};

static bool CompareHash (const unsigned char* digest, const char* hashStr) {
  char hashStrCmp[SHA_DIGEST_LENGTH * 2 + 1];
  for (int i = 0; i < SHA_DIGEST_LENGTH; ++i) {
    sprintf (&hashStrCmp[i * 2], "%.2x", digest[i]);
  }
  hashStrCmp[SHA_DIGEST_LENGTH * 2] = '\0';
  if (hashStr == hashStrCmp)
  {
    return true;
  }
  return false;
}

static void UpdateHashFromFrame (SFrameBSInfo& info, SHA1Context* ctx) {
  for (int i = 0; i < info.iLayerNum; ++i) {
    const SLayerBSInfo& layerInfo = info.sLayerInfo[i];
    int layerSize = 0;
    for (int j = 0; j < layerInfo.iNalCount; ++j) {
      layerSize += layerInfo.pNalLengthInByte[j];
    }
    SHA1Input (ctx, layerInfo.pBsBuf, layerSize);
  }
}

void* ThreadFunction(void* ptr)
{
  //------------------------------------------------------------
  // Get thread information
  igtl::MultiThreader::ThreadInfo* info = 
    static_cast<igtl::MultiThreader::ThreadInfo*>(ptr);

  //int id      = info->ThreadID;
  //int nThread = info->NumberOfThreads;
  ThreadData* td = static_cast<ThreadData*>(info->UserData);

  //------------------------------------------------------------
  // Get user data
  igtl::MutexLock::Pointer glock = td->glock;
  long interval = td->interval;
  std::cerr << "Interval = " << interval << " (ms)" << std::endl;
  //long interval = 1000;
  //long interval = (id + 1) * 100; // (ms)

  igtl::Socket::Pointer& socket = td->socket;

  //------------------------------------------------------------
  // Allocate TrackingData Message Class
  //
  // NOTE: TrackingDataElement class instances are allocated
  //       before the loop starts to avoid reallocation
  //       in each image transfer.

  ISVCEncoder* encoder_ = NULL;
  int rv = WelsCreateSVCEncoder (&encoder_);
  if (rv == 0 && encoder_ != NULL)
  {
    std::string fileName = polyFile;// + "/" + (std::string) kFileParamArray.pkcFileName;
    FileInputStream fileStream;
    fileStream.Open(fileName.c_str());
    SEncParamExt pEncParamExt;
    memset (&pEncParamExt, 0, sizeof (SEncParamExt));
    //kFileParamArray.iWidth =
    EncFileParamToParamExt (&kFileParamArray, &pEncParamExt);
    encoder_->InitializeExt(&pEncParamExt);
    int videoFormat = videoFormatI420;
    encoder_->SetOption (ENCODER_OPTION_DATAFORMAT, &videoFormat);
    while (!td->stop)
    {
      
      int frameSize = pEncParamExt.iPicWidth * pEncParamExt.iPicHeight * 3 / 2;
      
      BufferedData buf;
      buf.SetLength (frameSize);    
      SFrameBSInfo info;
      memset (&info, 0, sizeof (SFrameBSInfo));
      SSourcePicture pic;
      memset (&pic, 0, sizeof (SSourcePicture));
      SHA1Context ctx;
      memset (&ctx, 0, sizeof(SHA1Context));
      FILE* pFpBs = NULL;
      pFpBs = fopen ("encodedOutput", "wb");
      pic.iPicWidth    = pEncParamExt.iPicWidth;
      pic.iPicHeight   = pEncParamExt.iPicHeight;
      pic.iColorFormat = videoFormatI420;
      pic.iStride[0]   = pic.iPicWidth;
      pic.iStride[1]   = pic.iStride[2] = pic.iPicWidth >> 1;
      pic.pData[0]     = buf.data();
      pic.pData[1]     = pic.pData[0] + pEncParamExt.iPicWidth * pEncParamExt.iPicHeight;
      pic.pData[2]     = pic.pData[1] + (pEncParamExt.iPicWidth * pEncParamExt.iPicHeight >> 2);
      int iFrameIdx =0;
      while (fileStream.read (buf.data(), frameSize) == frameSize)
      {
        pic.uiTimeStamp = (long long)(iFrameIdx * (1000 / pEncParamExt.fMaxFrameRate));
        iFrameIdx++;
        //TestDebugCharArrayCmp(pic.pData[0], pic.pData[0], 200);
        int rv = encoder_->EncodeFrame (&pic, &info);
        //TestDebugCharArrayCmp(info.sLayerInfo[0].pBsBuf, info.sLayerInfo[0].pBsBuf, 200);
        if(rv == cmResultSuccess)
        {
          // 1. contain SHA encryption, could be removed, 2. contain the digest message could be as CRC
          UpdateHashFromFrame (info, &ctx);
          //---------------
          igtl::VideoMessage::Pointer videoMsg;
          videoMsg = igtl::VideoMessage::New();
          videoMsg->SetDeviceName("Video");
          videoMsg->SetBitStreamSize(info.iFrameSizeInBytes);
          videoMsg->AllocateScalars();
          videoMsg->SetScalarType(videoMsg->TYPE_UINT32);
          videoMsg->SetEndian(igtl_is_little_endian()==true?2:1); //little endian is 2 big endian is 1
          videoMsg->SetWidth(pic.iPicWidth);
          videoMsg->SetHeight(pic.iPicHeight);
          int frameSize = 0;
          int layerSize = 0;
          for (int i = 0; i < info.iLayerNum; ++i) {
            const SLayerBSInfo& layerInfo = info.sLayerInfo[i];
            layerSize = 0;
            for (int j = 0; j < layerInfo.iNalCount; ++j)
            {
              frameSize += layerInfo.pNalLengthInByte[j];
              layerSize += layerInfo.pNalLengthInByte[j];
            }
            //TestDebugCharArrayCmp(layerInfo.pBsBuf, layerInfo.pBsBuf, layerSize<200? layerSize:200);
            for (int i = 0; i < layerSize ; i++)
            {
              videoMsg->GetPackFragmentPointer(2)[frameSize-layerSize+i] = layerInfo.pBsBuf[i];
            }
            fwrite (layerInfo.pBsBuf, 1, layerSize, pFpBs); // write pure bit stream into file
          }
          //std::cerr<<"line break"<<std::endl;
          //TestDebugCharArrayCmp(videoMsg->GetPackFragmentPointer(2), videoMsg->GetPackFragmentPointer(1) + IGTL_VIDEO_HEADER_SIZE, 200);
          videoMsg->Pack();
          glock->Lock();
          
          for (int i = 0; i < videoMsg->GetNumberOfPackFragments(); i ++)
          {
            socket->Send(videoMsg->GetPackFragmentPointer(i), videoMsg->GetPackFragmentSize(i));
          }
          glock->Unlock();
          igtl::Sleep(interval);
        }
      }
      unsigned char digest[SHA_DIGEST_LENGTH];
      SHA1Result(&ctx, digest);
      CompareHash (digest, kFileParamArray.pkcHashStr);
      //------------------------------------------------------------
      // Loop
    }
  }
  WelsDestroySVCEncoder(encoder_);
  return NULL;
}


