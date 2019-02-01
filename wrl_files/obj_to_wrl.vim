function! Count( word )
    "https://vi.stackexchange.com/questions/6975/store-the-number-of-matches-in-vimscript-function
    redir => cnt
    silent exe '%s/' . a:word . '//gn'
    redir END

    let res = strpart(cnt, 0, stridx(cnt, " "))
    return res
endfunction
function! obj_to_wrl()
    "delete filler from blender, not generic for all obj's
    normal gg4dd
    " remove normals
    execute "g/vn /d"
    " delete random other junk
    execute "g/usemtl Material/d"
    execute "g/s off/d"
    
    let num_verts = Count("v ")
    let num_faces = Count("f ")
    execute "g/s off/d"
    " add vertex prefix
    normal gg
    execute "put!="
    normal gg
    let s = "points " . num_verts
    call setline('.', s)

    " add face prefix
    normal /f 
    execute "put!="
    let s2 = "coords " . num_faces
    call setline('.', s2)

    " remove v's and f's
    execute ":%s/v //g"
    execute ":%s/f //g"
    " remove normals per face
    " apparently " " and ' ' meann different things in vimscript
    execute ':%s/\(\d\+\)\/\/\(\d\+\)/\1/g'
    
    " reduce all obj vertex indices from base 1 to base 0
    normal gg
    execute 'g/\(\d\+\) \(\d\+\) \(\d\+\) \(\d\+\)'
    normal j
    normal VG
    "execute ":'<'>s/\(\zs\d\+\ze\)/\=(submatch(0)-1)/g"

    " add group ids
    normal G
    let s3 = "gids: " . num_faces
    execute "put="
    call setline('.', s3)
    execute "put="
    let s4 = repeat("0 ", num_faces)
    echo s4
    call setline('.', s4)
    execute "put="
    call setline('.', "groups 1")
    execute "put="
    call setline('.', "intpt 0 0 0 ")
    execute "put="
    call setline('.', "orient 1")


    " remove strange vim escape characters
    execute ':%s/\%x00//g'
endfunction
