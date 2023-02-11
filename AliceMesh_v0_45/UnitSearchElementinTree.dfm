object FormSearchElement_in_Tree: TFormSearchElement_in_Tree
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Search component in Tree'
  ClientHeight = 89
  ClientWidth = 593
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 0
    Top = 0
    Width = 222
    Height = 13
    Caption = 'Please, enter name element for search in Tree'
  end
  object EditNameElement_for_Search_in_Tree: TEdit
    Left = 0
    Top = 24
    Width = 593
    Height = 21
    TabOrder = 0
  end
  object ButtonOpenElement: TButton
    Left = 432
    Top = 64
    Width = 131
    Height = 25
    Caption = 'Open_Element'
    TabOrder = 1
    OnClick = ButtonOpenElementClick
  end
end
