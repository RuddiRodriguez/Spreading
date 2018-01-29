fid = fopen( 'myFile.txt', 'w' ) ;
 for cId = 1 : numel( vinterpto )
    fprintf( fid, '%f ', vinterpto{cId} ) ;
    fprintf( fid, '\n' ) ;
 end
 fclose( fid ) ;