unit Unitresidualplot;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, TeEngine, Series, ExtCtrls, TeeProcs, Chart, VclTee.TeeGDIPlus,
  Vcl.AppEvnts;

type
  TFormresidual = class(TForm)
    cht1: TChart;
    lnsrsSeries1: TLineSeries;
    lnsrsSeries2: TLineSeries;
    lnsrsSeries3: TLineSeries;
    lnsrsSeries4: TLineSeries;
    Timer1: TTimer;
    ApplicationEvents1: TApplicationEvents;
    procedure Timer1Timer(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
  private
    { Private declarations }
  public
    { Public declarations }
    brun_visible : Boolean;
  end;

var
  Formresidual: TFormresidual;

implementation

{$R *.dfm}

uses VisualUnit;

 // «апрет форме сворачиватьс€.
procedure TFormresidual.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

procedure TFormresidual.Timer1Timer(Sender: TObject);
var
   f : TStringList; // переменна€ типа объект TStringList
   i : Integer;
   fmin, fmax : Real;
   s, sub, subx : string;
begin
    // ƒействие будет происходить каждую секунду.
    f:=TStringList.Create();

      try
       if brun_visible then
       begin

       if (Laplas.egddata.itemper=0) then
     begin
     if (FileExists('statistic_convergence.txt')) then
      begin
         f.LoadFromFile('statistic_convergence.txt');

          if (FormatSettings.DecimalSeparator=',') then
          begin
             // заменить все точки в файле на зап€тые.
             for i:=0 to f.Count-1 do
             begin
                s:=f.Strings[i];
                f.Strings[i]:=StringReplace(s,'.',',',[rfReplaceAll]);
             end;
          end;

         // первые две строки нужно пропустить.
         Formresidual.cht1.SeriesList[0].Clear;
         Formresidual.cht1.SeriesList[1].Clear;
         Formresidual.cht1.SeriesList[2].Clear;
         Formresidual.cht1.SeriesList[3].Clear;
         for i:=2 to f.Count-1 do
         begin
            fmin:=20.0;
            fmax:=120.0;
            s:=Trim(f.Strings[i]);
            subx:=Trim(Copy(s,1,Pos(' ',s)));
            s:=Trim(Copy(s,Pos(' ',s),Length(s)));
            sub:=Trim(Copy(s,1,Pos(' ',s)));
            if (StrToFloat(sub)<fmin) then
            begin
               fmin:=StrToFloat(sub);
            end;
            if (StrToFloat(sub)>fmax) then
            begin
               fmax:=StrToFloat(sub);
            end;
            Formresidual.cht1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
            s:=Trim(Copy(s,Pos(' ',s),Length(s)));
            sub:=Trim(Copy(s,1,Pos(' ',s)));
            if (length(sub)>0) then
            begin
               if (StrToFloat(sub)<fmin) then
               begin
                  fmin:=StrToFloat(sub);
               end;
               if (StrToFloat(sub)>fmax) then
               begin
                  fmax:=StrToFloat(sub);
               end;
               Formresidual.cht1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
               s:=Trim(Copy(s,Pos(' ',s),Length(s)));
               sub:=Trim(Copy(s,1,Pos(' ',s)));
               if (length(sub)>0) then
               begin
                  if (StrToFloat(sub)<fmin) then
                  begin
                     fmin:=StrToFloat(sub);
                  end;
                  if (StrToFloat(sub)>fmax) then
                  begin
                     fmax:=StrToFloat(sub);
                  end;
                  Formresidual.cht1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clblue);
                  s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                  if (pos(' ',Trim(s))=0) then
                  begin
                    sub:=s;
                  end
                  else
                  begin
                     sub:=Trim(Copy(s,1,Pos(' ',s)));
                  end;
                  if (length(sub)>0) then
                  begin
                     if (StrToFloat(sub)<fmin) then
                     begin
                        fmin:=StrToFloat(sub);
                     end;
                     if (StrToFloat(sub)>fmax) then
                     begin
                        fmax:=StrToFloat(sub);
                     end;
                     Formresidual.cht1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
                  end
                  else
                  begin
                      // TODO
                      // обрыв данных после первых трЄх значений.
                  end;

                end
                else
                begin
                   // TODO
                   // обрыв данных после двух первых значений.
                end;
            end
            else
            begin
               // TODO
               // обрыв данных.
            end;
            Formresidual.cht1.LeftAxis.Minimum:=1e-3*fmin;
            Formresidual.cht1.LeftAxis.Maximum:=1e3*fmax;
         end;
         // Ќам ненужно запускать форму, нам нужно при запущенной
         // из вне формы посто€нно обновл€ть информацию.
         //Formresidual.Show;
      end
      else
      begin
         // ≈сли файл не найден на не надо ничего считывать посто€нно.
         //MainMemo.Lines.Add('statistic_convergence.txt not found.');
      end;
     end;
       end;

     except
        brun_visible:=false;

      end;

    f.Free;
end;

end.
