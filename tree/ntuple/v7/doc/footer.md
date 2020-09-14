# Layout

The main pieces are: 
* header, 
* footer, 
* and data.  

These are composed of various primitive and compound types. 

# Footer: 
Cluster metadata.  
```
* Frame 
* u64 (reserved) 
* u64 num clusters 
* Cluster summaries 
* Ntuple postscript 

footer-> +-+-+-+-+-+-+-+-+-+-+-+-+
         | Frame                 |
         +-+-+-+-+-+-+-+-+-+-+-+-+
         | (reserved)            |
         +-+-+-+-+-+-+-+-+-+-+-+-+
         | Num clusters          |
         +-+-+-+-+-+-+-+-+-+-+-+-+
         | Cluster summary [0]   |
         +-+-+-+-+-+-+-+-+-+-+-+-+
         |                       |
         .                       .
         +-+-+-+-+-+-+-+-+-+-+-+-+
         | Cluster summary [N-1] |
         +-+-+-+-+-+-+-+-+-+-+-+-+
         | Postscript            |
         +-+-+-+-+-+-+-+-+-+-+-+-+
```

## Postscript 
```
* u16 current version 
* u16 minimum version 
* u32 header size 
* u32 footer size 
* u32 CRC32 checksum 
=> 16 bytes 
```

# Primitive and compound types 

## Cluster 
```
* Uuid 
* Cluster summary 
* u32 num columns 
* [Column]
```

## ClusterSize 
```
* u32 num column elements
=> 4 bytes 
```

## Cluster summary 
```
* Frame 
* u64 cluster id 
* Version
* u64 first entry index 
* u64 num entries 
* Locator 
* u32 size
```

## Column 
``` 
* u64 column id 
* Column range 
* u32 num pages 
* [PageInfo]
```

## Column range
```
* u64 first element index 
* ClusterSize num column elements 
* i64 ROOT compression setting 
=> 20 bytes 
``` 

## Frame 
Version information.
```
* u16 current version 
* u16 minimum version 
* u32 frame size
=> 8 bytes 
``` 

## Locator 
```
i64 position
u32 bytes on storage 
String url
```

## PageInfo
```
* ClusterSize num column elements 
* Locator 
``` 

## String 
```
* u32 length 
* [u8] data 
```

## Uuid 
```
* Frame 
* String 
* u32 size
```

## Version
```
* Frame 
* u32 current version
* u32 minimum version
* u64 flags 
=> 16 bytes 
```
