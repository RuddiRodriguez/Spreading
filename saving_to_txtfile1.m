fid = fopen( 'myFile1.txt', 'w' ) ;
 for cId = 1 : numel( vinterpto )
     data = vinterpto{cId};
     data (data<0)=NaN;
    fprintf( fid, '%f ', data ) ;
    fprintf( fid, '\n' ) ;
 end
 fclose( fid ) ;