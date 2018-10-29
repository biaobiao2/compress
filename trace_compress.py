#encoding:utf-8

#===============================================================================
#    轨迹实时压缩、离线压缩算法。
#===============================================================================


import math
from datetime import datetime


EARTH_RADIUS = 6378137

def rad(pos):
    return pos * math.pi / 180.0;

#第一个点的经度、纬度，第二个点的经、度纬度
def getDistanceByPosition(alon,alat,blon,blat):
    
    alatrad = rad(alat)
    blatrad = rad(blat)
    
    diff_long = rad(alon) - rad(blon)
    e = math.acos( math.sin(blatrad)*math.sin(alatrad) + math.cos(alatrad)*math.cos(blatrad)*math.cos(diff_long) )
    
    e *= EARTH_RADIUS 
    return e

def getDistanceByTwoPosition(position1,position2):
    return getDistanceByPosition(position1[0],position1[1],position2[0],position2[1])


def point2line_distance(point,start_point,end_point):
    #任意一点到start与end直线的垂直距离
    
    start2end = getDistanceByTwoPosition(start_point,end_point)
    point2start = getDistanceByTwoPosition(point,start_point)
    point2end = getDistanceByTwoPosition(point,end_point)
    
    if start2end <= 10.0:
        #同一点
        return point2start
    
    p = (point2start + point2end + start2end)/2
    square = math.sqrt(p*(p-point2start)*(p-point2end)*(p-start2end))
    
    h = square*2/start2end
    
    return h

def top_down_compress_point(traces,start,end,threshold,trace_set):
    """ 采用 top_down算法迭代压缩轨迹：
    # 找出trace[start:end]中满足距离大于threshold的最大的一个点；若没有则返回None
    # 迭代求解
    
    Args:
        traces: 原始轨迹，由(time,long,lati)组成的list
        start：起始点序号
        end：结束点序号
        threshold：阈值
        trace_set：压缩后的轨迹点序号集合

    Returns:

    """
    if start+1 >= end:
        return
        
    point_distance_list = []
    for i in range(start+1,end):
        
        distance = point2line_distance(traces[i][1:],traces[start][1:],traces[end][1:])
        
        point_distance_list.append([i,distance])
    
    point_distance_list.sort(key=lambda d:d[1], reverse=True)
        
    if point_distance_list[0][1]>threshold:
        #print point_distance_list[0][0],point_distance_list[0][1]
        trace_set.add(point_distance_list[0][0])
        top_down_compress_point(traces,start,point_distance_list[0][0],threshold,trace_set)
        top_down_compress_point(traces,point_distance_list[0][0],end,threshold,trace_set)
    else:
        return        
    
def opening_window_compress_point(traces,threshold=10,remain_last_trace=False):
    '''实时轨迹压缩算法
    # 1.初始化一个容量为3、按时间有序的轨迹数据段；起点为第一个数据点、浮动点为第三个点；输出第一个点；
    # 2.如果每个点的误差（该点到由起点-浮动点组成的直线间垂直距离）都低于一个距离阈值，浮动点往后移动一位，转2；否则转3
    # 3.该点作为当前段的尾端点，下一段的起始点，浮动点往后移动一位；输出该点；转2
    # 4.计算复杂度 O(n^2)
    Args:
        traces: 原始轨迹
        threshold：距离阈值
        remain_last_trace：如果一段轨迹几乎不变，是否保留最后一个。原始opening window不保留，保留最后一次轨迹
    Returns:
        output_traces : 实时压缩后的轨迹.
    '''
    
    intial_traces = traces[0:3]
    
    #输出第一个点
    output_traces = []
    output_traces.append(intial_traces[0])
    print intial_traces[0]
    
    i = len(intial_traces)
    
    # 实时轨迹输入
    while i < len(traces)-1:
        
        shrift,new_traces,last_trace = should_shrift(intial_traces,threshold)
        
        if shrift:
            #输出

            if remain_last_trace is True and last_trace != output_traces[-1]:

                output_traces.append(last_trace)
                print last_trace
                
            output_traces.append(new_traces[0])
            print new_traces[0]
            
            new_traces.append(traces[i+1])
            intial_traces = new_traces
        else:
            intial_traces.append(traces[i+1])
            
        i += 1
        
    return output_traces

def should_shrift(traces,threshold):
    '''opening_window算法中判断是否traces中间的点到首尾直线的距离是否都低于阈值。
    # 若低于，返回(False,None)
    # 否则，返回新(True,new_traces)，new_traces的第一个点即为大于阈值的点
    '''
    for i in range(1,len(traces)-1):
        distance = point2line_distance(traces[i][1:],traces[0][1:],traces[-1][1:])
        
        if distance > threshold:
            return True,traces[i:],traces[i-1]
        
    return False,None,None
    
def compressTrace(traces):
    """
    Args:
        traces: 原始轨迹，由(time,long,lati)组成的list

    Returns:
        compress_traces : 压缩后的轨迹
    """

    compress_traces = []
    #记录起始点
    trace_set = set([0,len(traces)-1])
    #找出中间点
    top_down_compress_point(traces,0,len(traces)-1,300,trace_set)
     
    compressed_trace_index = sorted(trace_set)
    #print compressed_trace_index
    
    for index in compressed_trace_index:
        compress_traces.append(traces[index])
        #print traces[index][0],traces[index][1],traces[index][2]
        
    return compress_traces
        
def verifyTrace(traces):
    """
    Args:
        traces: 原始轨迹，由(time,long,lati)组成的list

    Returns:
        filtered_traces : 剔除有问题轨迹后的数据。判定依据：从后往前，比较相邻轨迹之间的时间差与距离
    """
    
    if len(traces) <= 1:
        return traces
    
    filtered_traces = [traces[-1],]
    
    for trace in traces[-2::-1]:
        #比较已存轨迹与历史轨迹之间的合理性。倒数第二个位置，逆序向前
        current_trace = filtered_traces[-1]
        
        #最大误差300m
        if abs(trace[1] - current_trace[1])<0.003 and abs(trace[2] - current_trace[2])<0.003:
            filtered_traces.append(trace)
        else:
            time_delta_seconds = (datetime.strptime(current_trace[0],"%Y-%m-%d_%H-%M-%S") - datetime.strptime(trace[0],"%Y-%m-%d_%H-%M-%S")).total_seconds()
            distance = getDistanceByTwoPosition((current_trace[1],current_trace[2]),(trace[1],trace[2]))
            # print trace,current_trace,distance,time_delta_seconds
            # 汽车理论最快速度 30m/s
            if distance/time_delta_seconds > 30:
                continue
            else:
                filtered_traces.append(trace)
            # 飞机、高铁等交通工具的位置跃变
            
    filtered_traces.reverse()
    
    return filtered_traces


